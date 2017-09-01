import numpy as np
import solid_state_tools as sst
import xml.etree.ElementTree as ET
import xml
from xml.dom import minidom
import subprocess
import periodictable as pt
import subprocess
import os
import time
import matplotlib.pyplot as plt
import re
import threading

plt.ion()

p_table = {i: el.__repr__() for i, el in enumerate(pt.elements)}
p_table_rev = {el.__repr__(): i for i, el in enumerate(pt.elements)}

hartree = 27.211


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
        self.engine_command = ["excitingser"]
        self.working_dirctory = '/exciting_files/'
        self.engine_process = None
        self.info_file = 'INFO.OUT'
        self.project_directory = None
        self.input_filename = 'input.xml'
        self.exciting_folder = self.find_exciting_folder()
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


        self.relax_file_timestamp = None


    def find_exciting_folder(self):
        p = subprocess.Popen(['which', 'excitingser'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        res, err = p.communicate()
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
            coords = np.array(self.split_and_remove_whitespace(basevec.text))
            crystal_base[i, :] = coords * scale

        atomic_cord_list = []
        for species in species_list:
            species_name = species.get('speciesfile').split('.')[0]
            species_number = p_table_rev[species_name]
            for j1, atom in enumerate(species):
                pos_vec = np.zeros(4)
                rel_coords = np.array(self.split_and_remove_whitespace(atom.get('coord')))
                if cart_bool:
                    pos_vec[:3] = np.dot(np.linalg.inv(crystal_base.T), rel_coords)
                else:
                    pos_vec[:3] = rel_coords
                pos_vec[3] = species_number
                atomic_cord_list.append(pos_vec)

        atom_array = np.array(atomic_cord_list)
        crystal_structure = sst.CrystalStructure(crystal_base, atom_array)
        return crystal_structure

    def make_tree(self):
        root = ET.Element("input")
        tree = ET.ElementTree(root)
        ET.SubElement(root, "title").text = self.general_options['title']
        return tree

    def add_scf_to_tree(self, tree, crystal_structure):
        root = tree.getroot()
        structure = ET.SubElement(root, "structure", speciespath=self.exciting_folder + 'species',tshift='false',autormt='true')
        crystal = ET.SubElement(structure, "crystal")
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

        groundstate = ET.SubElement(root, "groundstate", **self.scf_options)

    def add_bs_to_tree(self, tree, points):
        root = tree.getroot()
        properties = ET.SubElement(root, "properties")
        bandstructure = ET.SubElement(properties, "bandstructure")
        plot1d = ET.SubElement(bandstructure, "plot1d")
        path = ET.SubElement(plot1d, "path", **self.bs_options)
        for point in points:
            cords = point[0]
            label = point[1]
            ET.SubElement(path, "point", coord="{0:1.6f} {1:1.6f} {2:1.6f}".format(*cords), label=label)

    def add_relax_to_tree(self,tree):
        root = tree.getroot()
        relax = ET.SubElement(root, "relax",**self.relax_options)

    def add_gw_to_tree(self,tree,taskname='g0w0'):
        root = tree.getroot()
        relax = ET.SubElement(root, "gw",taskname=taskname,**self.gw_options)

    def start_gw(self,crystal_structure,band_structure_points=None):
        tree = self.make_tree()
        self.add_scf_to_tree(tree, crystal_structure)
        self.add_gw_to_tree(tree,taskname='g0w0')
        self.write_input_file(tree)

        def first_round():
            self.start_engine()
            self.engine_process.wait()
            if band_structure_points is not None:
                tree = self.make_tree()
                self.add_scf_to_tree(tree, crystal_structure)
                self.add_gw_to_tree(tree, taskname='band')
                self.add_bs_to_tree(tree,band_structure_points)
                self.write_input_file(tree)
                self.start_engine()

        t = threading.Thread(target=first_round)
        t.start()



    def start_relax(self,crystal_structure):
        tree = self.make_tree()
        self.add_scf_to_tree(tree, crystal_structure)
        self.add_relax_to_tree(tree)
        self.write_input_file(tree)
        time.sleep(0.05)
        self.start_engine()

    def load_relax_structure(self):
        file = self.project_directory+self.working_dirctory+'geometry_opt.xml'
        if not os.path.isfile(file):
            return None
        if self.relax_file_timestamp is not None and os.path.getmtime(file) == self.relax_file_timestamp:
            return None
        struc = self.parse_input_file(file)
        self.relax_file_timestamp = os.path.getmtime(file)
        return struc


    def write_input_file(self, tree):
        if not os.path.isdir(self.project_directory + self.working_dirctory):
            os.mkdir(self.project_directory + self.working_dirctory)
        xmlstr = minidom.parseString(ET.tostring(tree.getroot())).toprettyxml(indent="   ")
        with open(self.project_directory + self.working_dirctory + self.input_filename, "w") as f:
            f.write(xmlstr)

            # tree.write(self.project_directory+self.working_dirctory+self.input_filename)

    def start_engine(self):
        os.chdir(self.project_directory + self.working_dirctory)
        command = self.engine_command
        self.engine_process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        os.chdir(self.project_directory)

    def start_ground_state_calculation(self,crystal_structure,band_structure_points=None):
        tree = self.make_tree()
        self.add_scf_to_tree(tree, crystal_structure)
        if band_structure_points is not None:
            self.add_bs_to_tree(tree, band_structure_points)
        self.write_input_file(tree)
        time.sleep(0.05)
        self.start_engine()

    def kill_engine(self):
        if self.is_engine_running():
            self.engine_process.kill()

    def read_engine_status(self):
        if not self.is_engine_running():
            output, error = self.engine_process.communicate()
            return output, error
        else:
            return None

    def is_engine_running(self):
        if self.engine_process is None:
            return False
        if self.engine_process.poll() is None:
            return True
        else:
            return False

    def split_and_remove_whitespace(self, string):
        l1 = string.split()
        l2 = [float(x) for x in l1 if len(x) > 0]
        return l2

    def find_bandgap(self, bands):
        for i in range(len(bands)):
            valence_band = bands[i]
            cond_band = bands[i + 1]
            if (np.max(cond_band[:, 1]) > 0) and (np.min(cond_band[:, 1]) < 0):
                return None, None
            if any(cond_band[:, 1] > 0):
                break
        # Direct bandgap
        band_diff = cond_band[:, 1] - valence_band[:, 1]
        bandgap_index = np.argmin(band_diff)
        bandgap = np.min(band_diff)
        k_bandgap = valence_band[bandgap_index, 0]
        # Indirect bandgap
        if np.abs((np.min(cond_band[:, 1]) - np.max(valence_band[:, 1])) - bandgap) > 0.01:
            bandgap = (np.min(cond_band[:, 1]) - np.max(valence_band[:, 1]))
            k_bandgap = None
        return bandgap, k_bandgap

    def read_scf_status(self):
        try:
            f = open(self.project_directory + self.working_dirctory + self.info_file, 'r')
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
            e = xml.etree.ElementTree.parse(self.project_directory + self.working_dirctory + 'bandstructure.xml').getroot()
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

        bandgap, k_bandgap = self.find_bandgap(bands)
        special_k_points_together = zip(special_k_points, special_k_points_label)
        return sst.BandStructure(bands, bandgap, k_bandgap, special_k_points=special_k_points_together)

    def read_gw_bandstructure(self):
        band_structure = self.read_bandstructure()
        band_qp = np.loadtxt(self.project_directory+self.working_dirctory+'BAND-QP.OUT')
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

        bandgap,k_bandgap = self.find_bandgap(bands_qp)
        return sst.BandStructure(bands_qp,bandgap,k_bandgap,special_k_points=band_structure.special_k_points)


if __name__ == '__main__':
    os.chdir("/home/jannick/OpenDFT_projects/test/exciting_files")
    handler = Handler()
    handler.project_directory = "/home/jannick/OpenDFT_projects/test"

    print(handler.exciting_folder)
    tree = handler.make_tree()

    atoms = np.array([[0, 0, 0, 6], [0.25, 0.25, 0.25, 6]])
    unit_cell = 6.719 * np.array([[0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])

    crystal_structure = sst.CrystalStructure(unit_cell, atoms)
    handler.add_scf_to_tree(tree, crystal_structure)
    handler.add_bs_to_tree(tree, [[np.array([0, 0, 0]), "GAMMA"]])

    handler.write_input_file(tree)
    # plt.figure(2)
    #
    # handler.start_engine()
    # while handler.is_engine_running():
    #     time.sleep(1)
    #     res = handler.read_scf_status()
    #     if res is not None and len(res) >1:
    #         print('updating plot')
    #         plt.clf()
    #         plt.plot(res[:,0],res[:,1])
    #         plt.pause(0.01)
    #         plt.draw()
    # print('finished')
    # band_structure = handler.read_bandstructure()
    # plt.figure(1)
    # for band in band_structure.bands:
    #     plt.plot(band[:,0],band[:,1],color='b',linewidth=2)
    # for xc,xl in band_structure.special_k_points:
    #     plt.axvline(x=xc, color='k', linewidth=1.5)
    #
    # unzipped_k = zip(*band_structure.special_k_points)
    # plt.xticks(unzipped_k[0], unzipped_k[1], rotation='horizontal')
