import numpy as np
import solid_state_tools as sst
import xml.etree.ElementTree as ET
import xml
from xml.dom import minidom
import periodictable as pt
import subprocess
import os
import time
import re
import threading


p_table = {i: el.__repr__() for i, el in enumerate(pt.elements)}
p_table_rev = {el.__repr__(): i for i, el in enumerate(pt.elements)}

hartree = 27.211


if os.name == 'nt':
    shell_bool = True
else:
    shell_bool = False


def convert_greek(input):
    result = []
    for el in input:
        if el.lower() == 'gamma':
            result.append(r'$\Gamma$')
        else:
            result.append(el)
    return result


class Handler:
    def __init__(self):
        self.engine_name = 'exciting'
        self.default_extension = '.xml'
        self._engine_command = ["excitingser"]
        self._working_dirctory = '/exciting_files/'
        self.pseudo_dirctory = None
        self.engine_process = None
        self._info_file = 'INFO.OUT'
        self._filenames_tasks = {'scf': '/STATE.OUT', 'bandstructure': '/bandstructure.xml', 'g0w0 bands': '/BAND-QP.OUT',
                           'relax':'/geometry_opt.xml','g0w0':'/EIGVAL_GW.OUT','ks density':'/WF3D.xml'}

        self.supported_methods = sst.ComputationalMethods(['periodic', 'scf', 'gw', 'optical spectrum', 'phonons', 'relax'])

        self._timestamp_tasks = {}

        self.project_directory = None
        self._input_filename = 'input.xml'
        self.custom_command = ''
        self.custom_command_active = False
        self.dft_installation_folder = self.find_engine_folder()
        self.scf_options = {'do': 'fromscratch', 'nempty': '5', 'gmaxvr': '12.0', 'rgkmax': '7.0', 'ngridk': '1 1 1','frozencore':'false','xctype':'GGA_PBE'}
        self.scf_options_tooltip = {'do':r'Decides if the ground state is calculated starting from scratch, '
                                         'using the densities from file, or if its calculation is skipped and only the associated input parameters are read in.'}

        self.scf_options_tooltip['rgkmax'] = """The parameter rgkmax implicitly determines the number of basis functions and is one of the crucial parameters for the accuracy of the calculation. 
        It represents the product of two quantities: RMT,Min, the smallest of all muffin-tin radii, and |G+k|max, the maximum length for the G+k vectors. 
        Because each G+k vector represents one basis function, rgkmax gives the number of basis functions used for solving the Kohn-Sham equations. 
        Typical values of rgkmax are between 6 and 9. However, for systems with very short bond-lengths, significantly smaller values may be sufficient. 
        This may especially be the case for materials containing carbon, where rgkmax may be 4.5-5, or hydrogen, where even values between 3 and 4 may be sufficient. 
        In any case, a convergence check is indispensible for a proper choice of this parameter for your system!"""

        self.scf_options_tooltip['gmaxvr'] = r'Maximum length of G for expanding the interstitial density and potential.'
        self.scf_options_tooltip['xctype'] = """Type of exchange-correlation functional to be used

    No exchange-correlation funtional ( Exc=0 )
    LDA, Perdew-Zunger/Ceperley-Alder, Phys. Rev. B 23, 5048 (1981)
    LSDA, Perdew-Wang/Ceperley-Alder, Phys. Rev. B 45, 13244 (1992)
    LDA, X-alpha approximation, J. C. Slater, Phys. Rev. 81, 385 (1951)
    LSDA, von Barth-Hedin, J. Phys. C 5, 1629 (1972)
    GGA, Perdew-Burke-Ernzerhof (PBE), Phys. Rev. Lett. 77, 3865 (1996)
    GGA, Revised PBE, Zhang-Yang, Phys. Rev. Lett. 80, 890 (1998)
    GGA, PBEsol, arXiv:0707.2088v1 (2007)
    GGA, asymptotically corrected PBE (acPBE), arXiv:1409.4834 (2014)
    GGA, Wu-Cohen exchange (WC06) with PBE correlation, Phys. Rev. B 73, 235116 (2006)
    GGA, Armiento-Mattsson (AM05) spin-unpolarised functional, Phys. Rev. B 72, 085108 (2005)
    EXX, Exact Exchange, Phys. Rev. Lett. 95, 136402 (2005)
    Hybrid, PBE0, J. Chem. Phys. 110, 5029 (1999)

Type: 	choose from:
LDA_PZ
LDA_PW
LDA_XALPHA
LDA_vBH
GGA_PBE
GGA_PBE_R
GGA_PBE_SOL
GGA_WC
GGA_AM05
GGA_AC_PBE
HYB_PBE0
HYB_LDA0
EXX
none
Default: 	GGA_PBE"""

        self.scf_options_tooltip['ngridk'] = """Number of k grid points along the basis vector directions. 
        Alternatively give autokpt and radkpt, or nktot. In the latter cases any value given for ngridk is not used. 
        
        Notes: Phonon calculations using supercells adjust the k-grid according to the supercell size; if the element xs is given, 
        the present attribute is overwritten by the value in xs for xs-related groundstate calculations; 
        the values of the present attribute are also relevant for calculations related to the element gw."""

        self.scf_options_tooltip['nempty'] = """Defines the number of eigenstates beyond that required for charge neutrality. 
        When running metals it is not known a priori how many states will be below the Fermi energy for each k-point. 
        Setting nempty greater than zero allows the additional states to act as a buffer in such cases. 
        Furthermore, magnetic calculations use the first-variational eigenstates as a basis for setting up the second-variational Hamiltonian, 
        and thus nempty will determine the size of this basis set. Convergence with respect to this quantity should be checked."""

        self.general_options = {'title': 'title'}
        self.bs_options = {'steps': '300'}
        self.relax_options = {'method':'bfgs'}

        self.gw_options = {'nempty':'0','ngridq':"2 2 2",'ibgw':'1','nbgw':'0'}

        self.gw_options_tooltip = {'nempty':'Number of empty states (cutoff parameter) used in GW.\n'
                                            'If not specified, the same number as for the groundstate calculations is used.'}
        self.gw_options_tooltip['ngridq'] = 'k/q-point grid size to be used in GW calculations'
        self.gw_options_tooltip['ibgw'] = 'Lower band index for GW output.'
        self.gw_options_tooltip['nbgw'] = 'Upper band index for GW output.'

        self.phonons_options = {'do':'fromscratch','ngridq':'2 2 2'}

        self.optical_spectrum_options = {'xstype':'BSE','intv':'-0.5 0.5','points':'1000','bsetype':"singlet",'nstlbse':"1 4 1 4",'screentype':'full'
                                         ,'nempty_screeing':'0','use gw':'false','nempty':'5','ngridq':"4 4 4",'ngridk':"4 4 4",'broad':'0.01',
                                         'gqmax':"0.0",'vkloff':"0.0 0.0 0.0"}

        self.relax_file_timestamp = None

    def find_engine_folder(self):
        p = subprocess.Popen(['which', 'excitingser'], stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=shell_bool)
        res, err = p.communicate()
        res = res.decode()
        res = res.split('bin')[0]
        return res.strip()

    def parse_input_file(self, filename):
        tree = ET.parse(filename)
        root = tree.getroot()

        structure = root.find('structure')
        crystal = structure.find('crystal')
        species_list = structure.findall('species')

        cart_bool = structure.get('cartesian')
        if cart_bool is None:
            cart_bool = False
        else:
            cart_bool = bool(cart_bool)

        basevecs = crystal.findall('basevect')
        try:
            scale = float(crystal.get('scale'))
        except Exception:
            scale = 1
        crystal_base = np.zeros((3, 3))

        for i, basevec in enumerate(basevecs):
            coords = np.array(self._split_and_remove_whitespace(basevec.text))
            crystal_base[i, :] = coords * scale

        atomic_cord_list = []
        for species in species_list:
            species_name = species.get('speciesfile').split('.')[0]
            species_number = p_table_rev[species_name]
            for j1, atom in enumerate(species):
                pos_vec = np.zeros(4)
                rel_coords = np.array(self._split_and_remove_whitespace(atom.get('coord')))
                if cart_bool:
                    pos_vec[:3] = np.dot(np.linalg.inv(crystal_base.T), rel_coords)
                else:
                    pos_vec[:3] = rel_coords
                pos_vec[3] = species_number
                atomic_cord_list.append(pos_vec)

        atom_array = np.array(atomic_cord_list)
        crystal_structure = sst.CrystalStructure(crystal_base, atom_array,scale=scale)
        return crystal_structure

    def start_ground_state(self, crystal_structure, band_structure_points=None):
        try:
            os.remove(self.project_directory + self._working_dirctory + '/INFO.OUT')
        except Exception as e:
            print(e)
        self._read_timestamps()
        tree = self._make_tree()
        self._add_scf_to_tree(tree, crystal_structure)
        if band_structure_points is not None:
            self._add_bs_to_tree(tree, band_structure_points)
        self._write_input_file(tree)
        time.sleep(0.05)
        self._start_engine()

    def start_optical_spectrum(self,crystal_structure):
        self._filenames_tasks['optical spectrum'] = '/EPSILON_BSE' + self.optical_spectrum_options['bsetype'] + '_SCRfull_OC11.OUT'
        self._read_timestamps()
        tree = self._make_tree()
        self._add_scf_to_tree(tree, crystal_structure)
        self._add_optical_spectrum_to_tree(tree)
        self._write_input_file(tree)
        time.sleep(0.05)
        self._start_engine()

    def start_gw(self,crystal_structure,band_structure_points=None):
        self._read_timestamps()
        tree = self._make_tree()
        self._add_scf_to_tree(tree, crystal_structure)
        self._add_gw_to_tree(tree, taskname='g0w0')
        self._write_input_file(tree)

        def first_round():
            self._start_engine()
            if self.custom_command_active:
                while self.is_engine_running(tasks=['g0w0']):
                    time.sleep(1)
            else:
                self.engine_process.wait()
            if band_structure_points is not None:
                tree = self._make_tree()
                self._add_scf_to_tree(tree, crystal_structure)
                self._add_gw_to_tree(tree, taskname='band')
                self._add_bs_to_tree(tree, band_structure_points)
                self._write_input_file(tree)
                self._start_engine()

        t = threading.Thread(target=first_round)
        t.start()

    def start_phonon(self, crystal_structure, band_structure_points):
        self._read_timestamps()
        tree = self._make_tree()
        self._add_scf_to_tree(tree, crystal_structure)
        self._add_phonon_to_tree(tree, band_structure_points)
        self._write_input_file(tree)
        time.sleep(0.05)
        self._start_engine()

    def start_relax(self,crystal_structure):
        self._read_timestamps()
        tree = self._make_tree()
        self._add_scf_to_tree(tree, crystal_structure)
        self._add_relax_to_tree(tree)
        self._write_input_file(tree)
        time.sleep(0.05)
        self._start_engine()

    def load_relax_structure(self):
        file = self.project_directory+self._working_dirctory + 'geometry_opt.xml'
        if not os.path.isfile(file):
            return None
        if self.relax_file_timestamp is not None and os.path.getmtime(file) == self.relax_file_timestamp:
            return None
        struc = self.parse_input_file(file)
        self.relax_file_timestamp = os.path.getmtime(file)
        return struc

    def read_scf_status(self):
        try:
            f = open(self.project_directory + self._working_dirctory + self._info_file, 'r')
        except IOError:
            return None
        info_text = f.read()
        f.close()
        scf_list = []
        matches = re.findall(r"SCF iteration number[\s\t]*:[\s\t]*\d*", info_text)
        for match in matches:
            ms = match.split(':')
            scf_list.append(int(ms[1]))

        scf_energy_list = []
        matches = re.findall(r"Total energy [\s\t]*:[\s\t]*[-+]?\d*\.\d+", info_text)
        for match in matches:
            ms = match.split(':')
            scf_energy_list.append(float(ms[1]))

        res = np.array(zip(scf_list, scf_energy_list))
        if len(res) < 2:
            return None
        return res

    def read_bandstructure(self):
        try:
            e = xml.etree.ElementTree.parse(self.project_directory + self._working_dirctory + 'bandstructure.xml').getroot()
        except IOError as e:
            return None

        titlestring = e.find('title').itertext().next()

        bands = []

        for atype in e.findall('band'):
            band = []
            for point in atype.findall('point'):
                band.append([point.get('distance'), point.get('eval')])
            bands.append(np.array(band, dtype=np.float))

        special_k_points = []
        special_k_points_label = []
        for vertex in e.findall('vertex'):
            special_k_points.append(float(vertex.get('distance')))
            special_k_points_label.append(vertex.get('label'))

        special_k_points_label = convert_greek(special_k_points_label)

        empirical_correction = 0

        for band in bands:
            band[:, 1] = band[:, 1] * hartree
            if empirical_correction != 0:
                if (band[:, 1] > 0).any():
                    band[:, 1] = band[:, 1] + empirical_correction / 2
                else:
                    band[:, 1] = band[:, 1] - empirical_correction / 2

        special_k_points_together = zip(special_k_points, special_k_points_label)


        return sst.BandStructure(bands, special_k_points=special_k_points_together)

    def read_gw_bandstructure(self,filename='BAND-QP.OUT'):
        band_structure = self.read_bandstructure()
        band_qp = np.loadtxt(self.project_directory + self._working_dirctory + filename)
        k_qp = band_qp[:, 0]
        E_qp = band_qp[:, 1]*hartree

        band_sep = np.where(k_qp == 0)[0]
        N_k = band_sep[1]
        N_bands = len(band_sep)
        k_qp = np.resize(k_qp, [N_bands, N_k]).T
        E_qp = np.resize(E_qp, [N_bands, N_k]).T
        bands_qp = []
        for i in range(N_bands):
            bands_qp.append(np.array([k_qp[:,i],E_qp[:,i]]).T )

        bandgap,k_bandgap = self._find_bandgap(bands_qp)
        if filename == 'PHDISP.OUT':
            bs_type = 'phonon'
        else:
            bs_type = 'gw'
        return sst.BandStructure(bands_qp,bandgap,k_bandgap,special_k_points=band_structure.special_k_points)

    def read_phonon_bandstructure(self):

        return self.read_gw_bandstructure(filename='PHDISP.OUT')

    def read_optical_spectrum(self):
        eps_11 = np.loadtxt(self.project_directory + self._working_dirctory + 'EPSILON_BSE' + self.optical_spectrum_options['bsetype'] + '_SCRfull_OC11.OUT')
        eps_22 = np.loadtxt(self.project_directory + self._working_dirctory + 'EPSILON_BSE' + self.optical_spectrum_options['bsetype'] + '_SCRfull_OC22.OUT')
        eps_33 = np.loadtxt(self.project_directory + self._working_dirctory + 'EPSILON_BSE' + self.optical_spectrum_options['bsetype'] + '_SCRfull_OC33.OUT')


        list_of_eps2 = [eps_11[:,2],eps_22[:,2],eps_33[:,2]]
        list_of_eps1 = [eps_11[:,1],eps_22[:,1],eps_33[:,1]]


        return sst.OpticalSpectrum(eps_11[:,0]*hartree,list_of_eps2,epsilon1=list_of_eps1)

    def read_ks_state(self):
        self._convert_3d_plot()
        l_data = np.genfromtxt(self.project_directory + self._working_dirctory + 'WF3D.xsf', skip_header=9, skip_footer=2, dtype=np.float)
        n_grid = l_data.shape[1]
        # data = l_data.reshape((l_data.shape[1],l_data.shape[1],l_data.shape[1]),order='C')
        data = np.zeros((n_grid,n_grid,n_grid))
        for i in range(n_grid):
            for j in range(n_grid):
                for k in range(n_grid):
                    data[i,k,j] = l_data[n_grid*j+k,i]
            # data[:,:,i] = l_data[:,i].reshape((n_grid,n_grid),order='F')
        data = data/data.max()
        return sst.KohnShamDensity(data)

    def calculate_ks_density(self,crystal_structure,bs_point,grid='40 40 40'):
        self._read_timestamps()
        tree = self._make_tree()
        self._add_scf_to_tree(tree, crystal_structure, skip=True)
        self._add_ks_density_to_tree(tree, bs_point, grid)
        self._write_input_file(tree)
        self._start_engine()

    def calculate_electron_density(self, crystal_structure):
        # TODO implement this
        raise NotImplementedError

    def kill_engine(self):
        try:
            self.engine_process.kill()
        except Exception as e:
            print(e)

    def is_engine_running(self,tasks=None):
        if self.custom_command_active:
            return self._is_engine_running_custom_command(tasks)
        else:
            if self.engine_process is None:
                return False
            if self.engine_process.poll() is None:
                return True
            else:
                return False

    def will_scf_run(self):
        if self.scf_options['do'] == 'skip':
            return False
        else:
            return True

    def _make_tree(self):
        root = ET.Element("input")
        tree = ET.ElementTree(root)
        ET.SubElement(root, "title").text = self.general_options['title']
        return tree

    def _add_scf_to_tree(self, tree, crystal_structure, skip=False):
        root = tree.getroot()
        structure = ET.SubElement(root, "structure", speciespath=self.dft_installation_folder + 'species', tshift='false', autormt='true')
        crystal = ET.SubElement(structure, "crystal",scale='1.0')
        for i in range(3):
            lattice_vector = crystal_structure.lattice_vectors[i, :]
            ET.SubElement(crystal, "basevect").text = "{0:1.6f} {1:1.6f} {2:1.6f}".format(*lattice_vector)

        abs_coord_atoms = crystal_structure.atoms
        species = set(abs_coord_atoms[:, 3].astype(np.int))
        n_species = len(species)
        n_atoms = abs_coord_atoms.shape[0]

        for specie in species:
            species_mask = abs_coord_atoms[:, 3].astype(np.int) == specie
            sub_coords = abs_coord_atoms[species_mask, :]
            specie_xml_el = ET.SubElement(structure, "species", speciesfile=p_table[specie] + '.xml')
            n_atoms_specie = sub_coords.shape[0]
            for i in range(n_atoms_specie):
                ET.SubElement(specie_xml_el, "atom", coord="{0:1.6f} {1:1.6f} {2:1.6f}".format(*sub_coords[i, :]))


        if skip:
            new_scf_options = {}
            new_scf_options.update(self.scf_options)
            new_scf_options['do'] = 'skip'
            groundstate = ET.SubElement(root, "groundstate",**new_scf_options)
        else:
            groundstate = ET.SubElement(root, "groundstate",**self.scf_options)

    def _add_bs_to_tree(self, tree, points):
        root = tree.getroot()
        properties = ET.SubElement(root, "properties")
        bandstructure = ET.SubElement(properties, "bandstructure")
        plot1d = ET.SubElement(bandstructure, "plot1d")
        path = ET.SubElement(plot1d, "path", **self.bs_options)
        for point in points:
            cords = point[0]
            label = point[1]
            ET.SubElement(path, "point", coord="{0:1.6f} {1:1.6f} {2:1.6f}".format(*cords), label=label)

    def _add_relax_to_tree(self, tree):
        root = tree.getroot()
        relax = ET.SubElement(root, "relax",**self.relax_options)

    def _add_gw_to_tree(self, tree, taskname='g0w0'):
        root = tree.getroot()
        gw = ET.SubElement(root, "gw",taskname=taskname,**self.gw_options)

    def _add_optical_spectrum_to_tree(self, tree):
        root = tree.getroot()
        if self.optical_spectrum_options['use gw'].lower() == 'true':
            gw = ET.SubElement(root, "gw", taskname='skip')
        elif self.optical_spectrum_options['use gw'].lower() == 'false':
            pass
        else:
            raise Exception('Bad option for use gw')

        xs = ET.SubElement(root, "xs", xstype = self.optical_spectrum_options['xstype'],nempty = self.optical_spectrum_options['nempty'],
                           ngridk=self.optical_spectrum_options['ngridk'],ngridq = self.optical_spectrum_options['ngridq']
                           ,broad=self.optical_spectrum_options['broad'],gqmax=self.optical_spectrum_options['gqmax'],
                           vkloff=self.optical_spectrum_options['vkloff'])
        ewindow = ET.SubElement(xs, "energywindow",intv = self.optical_spectrum_options['intv'],points = self.optical_spectrum_options['points'])
        screening = ET.SubElement(xs, "screening",screentype = self.optical_spectrum_options['screentype'],
                                  nempty=self.optical_spectrum_options['nempty_screeing'])
        BSE = ET.SubElement(xs, "BSE", bsetype = self.optical_spectrum_options['bsetype'],nstlbse = self.optical_spectrum_options['nstlbse'])
        qpointset = ET.SubElement(xs,'qpointset')
        qpoint = ET.SubElement(qpointset,'qpoint').text = '0.0 0.0 0.0'

    def _add_phonon_to_tree(self, tree, points):
        root = tree.getroot()
        phonon = ET.SubElement(root, "phonons",**self.phonons_options)
        disp = ET.SubElement(phonon, "phonondispplot")
        plot1d = ET.SubElement(disp, "plot1d")
        path = ET.SubElement(plot1d,'path', steps=self.bs_options['steps'])
        for point in points:
            cords = point[0]
            label = point[1]
            ET.SubElement(path, "point", coord="{0:1.6f} {1:1.6f} {2:1.6f}".format(*cords), label=label)

    def _read_timestamps(self):
        for task,filename in self._filenames_tasks.items():
            if os.path.isfile(self.project_directory+self._working_dirctory+filename):
                self._timestamp_tasks[task]=os.path.getmtime(self.project_directory + self._working_dirctory + filename)

    def _check_if_scf_is_finished(self):
        try:
            f = open(self.project_directory + self._working_dirctory + self._info_file, 'r')
        except IOError:
            return False
        info_text = f.read()
        f.close()
        scf_list = []
        matches = re.findall(r"EXCITING[\s\w]*stopped", info_text)
        if len(matches) == 0:
            return False
        else:
            return True

    def _write_input_file(self, tree):
        if not os.path.isdir(self.project_directory + self._working_dirctory):
            os.mkdir(self.project_directory + self._working_dirctory)
        xmlstr = minidom.parseString(ET.tostring(tree.getroot())).toprettyxml(indent="   ")
        with open(self.project_directory + self._working_dirctory + self._input_filename, "w") as f:
            f.write(xmlstr)

            # tree.write(self.project_directory+self.working_dirctory+self.input_filename)

    def _start_engine(self):
        os.chdir(self.project_directory + self._working_dirctory)
        if self.custom_command_active:
            command = ['bash',self.custom_command]
            # if tasks is not None:
            #     filenames = [self.filenames_tasks[task] for task in tasks if task != 'scf']
            #     for filename in filenames:
            #         os.remove(self.project_directory+self.working_dirctory+filename)
        else:
            command = self._engine_command

        self.engine_process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        os.chdir(self.project_directory)

    #
    # def read_engine_status(self):
    #     if not self.is_engine_running():
    #         output, error = self.engine_process.communicate()
    #         return output, error
    #     else:
    #         return None

    def _is_engine_running_custom_command(self,tasks):
        bool_list = []

        for task in tasks:
            if task == 'scf':
                continue
            filename = self._filenames_tasks[task]
            # filenames[task] = filename
            file_exists = os.path.isfile(self.project_directory + self._working_dirctory + filename)
            # files_exists[task] = file_exists
            if file_exists:
                timestamp = os.path.getmtime(self.project_directory + self._working_dirctory + filename)
                try:
                    file_is_old_bool = timestamp == self._timestamp_tasks[task]
                except KeyError:
                    file_is_old_bool = False

            file_exist_and_new = file_exists and not file_is_old_bool
            bool_list.append(file_exist_and_new)


        if 'scf' in tasks:
            bool_list.append(self._check_if_scf_is_finished())

        finished_bool = all(bool_list)
        if finished_bool:
            return False
        else:
            return True

    def _split_and_remove_whitespace(self, string):
        l1 = string.split()
        l2 = [float(x) for x in l1 if len(x) > 0]
        return l2

    def _convert_3d_plot(self):
        os.chdir(self.project_directory + self._working_dirctory)
        command = 'xsltproc $EXCITINGVISUAL/plot3d2xsf.xsl WF3D.xml > WF3D.xsf'
        self.helper_process = subprocess.call(command,shell=True)

    def _add_ks_density_to_tree(self, tree, bs_point, grid):
        root = tree.getroot()
        properties = ET.SubElement(root, "properties")
        wfplot = ET.SubElement(properties, "wfplot")
        kstlist = ET.SubElement(wfplot, "kstlist")
        pointstatepair = ET.SubElement(kstlist, "pointstatepair").text = '{0} {1}'.format(*bs_point)

        plot3d = ET.SubElement(wfplot, "plot3d")
        box = ET.SubElement(plot3d, "box",grid=grid)

        p0 = np.array([0.0,0.0,0.0])
        p1 = np.array([1.0,0.0,0.0])
        p2 = np.array([0.0,1.0,0.0])
        p3 = np.array([0.0,0.0,1.0])


        origin = ET.SubElement(box, "origin",coord="{0:1.6f} {1:1.6f} {2:1.6f}".format(*p0))
        p1 = ET.SubElement(box, "point",coord="{0:1.6f} {1:1.6f} {2:1.6f}".format(*p1))
        p2 = ET.SubElement(box, "point",coord="{0:1.6f} {1:1.6f} {2:1.6f}".format(*p2))
        p3 = ET.SubElement(box, "point",coord="{0:1.6f} {1:1.6f} {2:1.6f}".format(*p3))


if __name__ == '__main__':
    handler = Handler()
    handler.project_directory = "/home/jannick/OpenDFT_projects/visualize"

    handler.custom_command_active = True
    ac_bo  = handler.is_engine_running(tasks = ['scf','bandstructure','g0w0'])
    ac_bo2 = handler._check_if_scf_is_finished()

    print(ac_bo,ac_bo2)

    # print(handler.exciting_folder)
    # tree = handler.make_tree()
    #
    # atoms = np.array([[0, 0, 0, 6], [0.25, 0.25, 0.25, 6]])
    # unit_cell = 6.719 * np.array([[0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    #
    # crystal_structure = sst.CrystalStructure(unit_cell, atoms)
    # handler.calculate_ks_density(crystal_structure,[1,1])
    # KS_dens = handler.load_ks_state()
    #
    #
    # from mayavi import mlab
    #
    # unit_cell = 6.719 * np.array([[0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    #
    # for i1, a in enumerate(unit_cell):
    #     i2 = (i1 + 1) % 3
    #     i3 = (i1 + 2) % 3
    #     for b in [np.zeros(3), unit_cell[i2]]:
    #         for c in [np.zeros(3), unit_cell[i3]]:
    #             p1 = b + c
    #             p2 = p1 + a
    #             mlab.plot3d([p1[0], p2[0]],
    #                         [p1[1], p2[1]],
    #                         [p1[2], p2[2]],
    #                         tube_radius=0.1)
    #
    # cp = mlab.contour3d(KS_dens.density, contours=10, transparent=True,
    #                     opacity=0.5, colormap='hot')
    # # Do some tvtk magic in order to allow for non-orthogonal unit cells:
    # polydata = cp.actor.actors[0].mapper.input
    # pts = np.array(polydata.points) - 1
    # # Transform the points to the unit cell:
    # polydata.points = np.dot(pts, unit_cell / np.array(KS_dens.density.shape)[:, np.newaxis])
    #
    #
    # mlab.view(azimuth=155, elevation=70, distance='auto')
    #
    # mlab.show()


