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
    search_command = 'where'
else:
    shell_bool = False
    search_command = 'which'


def convert_greek(input):
    result = []
    for el in input:
        if el.lower() == 'gamma':
            result.append(r'$\Gamma$')
        else:
            result.append(el)
    return result


class Handler:
    """Main class for interacting with the electronic structure calculation engine EXCITING."""
    def __init__(self):
        self.engine_name = 'exciting'
        self.default_extension = '.xml'
        self._engine_command = ["excitingser"]
        self.working_dirctory = '/exciting_files/'
        self.pseudo_directory = None
        self.engine_process = None
        self.info_file = 'INFO.OUT'
        self.info_text = """<a href="http://exciting-code.org/">exciting</a> is an all-electron full-potential computer package <a href="http://iopscience.iop.org/0953-8984/26/36/363202">[GUL-2014]</a> for first-principles calculations, based on (linearized) augmented planewave + local orbital [(L)APW+lo] methods. 
        This family of basis sets is known as the most precise numerical scheme to solve the Kohn-Sham equations of density-functional theory (DFT), reaching extremely high - up to muHartree - precision <a href="http://iopscience.iop.org/0953-8984/26/36/363202">[GUL-2014]</a>, <a href="http://science.sciencemag.org/content/sci/351/6280/aad3000.full.pdf?ijkey=teUZMpwU49vhY&keytype=ref&siteid=sci">[LEJ-2016]</a>. Different schemes are available to account for van der Waals forces.

As suggested by its name, exciting has a major focus on excited-state properties. It includes a module for time-dependent DFT (TDDFT) in the linear-response regime [SAG-2009], implementing a number of static and dynamical exchange-correlation kernels[SAG-2009], [RIG-2015]. TDDFT is preferably adopted to compute absorption and electron-loss spectra of materials with weak electron-hole interaction, such as small molecules and metals, also at finite momentum transfer [ALK-2013]. 
For systems with pronounced correlation effects, exciting offers a rich spectroscopy module based on many-body perturbation theory. 
To compute quasi-particle band structures, the GW approach is implemented in the single-shot G0W0 approximation [NAB-2016]. 
Recent developments of the PBE0 hybrid functional and of the LDA-1/2 method provide improved starting points for G0W0 calculations [PEL-2016]. 
The solution of the Bethe-Salpeter equation (BSE) offers an accurate description of excitations in the valence [SAG-2009] and in the core region [VOR-2017] on the same footing. 
Specific modules of exciting are dedicated to advanced light-matter interaction processes, such as Raman scattering, second-harmonic generation [SHA-2004], and the magneto-optic Kerr effect [GUL-2014].

exciting is an open-source, GPL-licensed code, written in a clean and fully documented programming style, with a modern source-code management, a dynamical build system, and automated tests. 
It is equipped with a detailed documentation of current developments, including an interactive Input Reference webpage and over 30 Tutorials illustrating basic and advanced features. 
The interface with pre- and post-processing tools integrates the capabilities of exciting for specific tasks, such as calculating elastic constants [GOL-2013] and optical coefficients [VOR-2016], as well as performing a cluster expansion for, e.g., thermoelectric materials with large parent cells [TRO-2017]."""

        self._filenames_tasks = {'scf': '/STATE.OUT', 'bandstructure': '/bandstructure.xml', 'g0w0 bands': '/BAND-QP.OUT',
                           'relax':'/geometry_opt.xml','g0w0':'/EIGVAL_GW.OUT','ks density':'/WF3D.xml','optical spectrum': None}

        self.supported_methods = sst.ComputationalMethods(['periodic', 'scf', 'g0w0', 'optical spectrum', 'phonons', 'relax','bandstructure'])

        self._timestamp_tasks = {}

        self.project_directory = None
        self.input_filename = 'input.xml'

        self.current_input_file = self.input_filename
        self.current_output_file = self.info_file

        self.custom_command = ''
        self.custom_command_active = False
        self.dft_installation_folder = self.find_engine_folder()
        self.scf_options = {'do': 'fromscratch', 'nempty': '5', 'gmaxvr': '12.0', 'rgkmax': '5.0', 'ngridk': '5 5 5','frozencore':'false','xctype':'GGA_PBE'}
        self.scf_options_tooltip = {'do':"""Decides if the ground state is calculated starting from scratch, using the densities from file, or if its calculation is skipped and only the associated input parameters are read in.
Type: 	choose from:
fromscratch
fromfile
skip"""}

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
        This parameter is critical for the convergence of the results!
        Alternatively give autokpt and radkpt, or nktot. In the latter cases any value given for ngridk is not used. 
        
        For molecules you can use 1 1 1. For medium size cells something around 4 4 4 should suffice. For very small unit cells, like e.g. diamond, something like 7 7 7
        should be ok.  
        
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
        self.relax_options_tooltip = {}

        self.gw_options = {'nempty':'0','ngridq':"2 2 2",'ibgw':'1','nbgw':'0'}

        self.gw_options_tooltip = {'nempty':'Number of empty states (cutoff parameter) used in GW.\n'
                                            'If not specified, the same number as for the groundstate calculations is used.'}
        self.gw_options_tooltip['ngridq'] = 'k/q-point grid size to be used in GW calculations'
        self.gw_options_tooltip['ibgw'] = 'Lower band index for GW output.'
        self.gw_options_tooltip['nbgw'] = 'Upper band index for GW output.'

        self.phonons_options = {'do':'fromscratch','ngridq':'2 2 2'}
        self.phonons_options_tooltip = {}

        self.optical_spectrum_options = {'xstype':'BSE','intv':'0.0 0.5','points':'1000','bsetype':"singlet",'nstlbse':"1 4 1 4",'screentype':'full'
                                         ,'nempty_screening':'0','use gw':'false','nempty':'5','ngridq':"4 4 4",'ngridk':"4 4 4",'broad':'0.0005',
                                         'gqmax':"3.0",'vkloff':"0.231 0.151 0.432"}
        self.optical_spectrum_options_tooltip = {'nstlbse':'Range of bands included for the BSE calculation.\nThe first pair of numbers corresponds to the band index for local orbitals and valence states (counted from the lowest eigenenergy),\nthe second pair corresponds to the band index of the conduction states (counted from the Fermi level).',
                                                 'nempty_screening':'Number of empty states.',
                                                 'vkloff':'k-point offset for screening Type: floats\n Do not use negative numbers',
                                                 'intv':'energy interval lower and upper limits (in hartree)',
                                                 'gqmax':'|G+q| cutoff for Kohn-Sham response function, screening and for expansion of Coulomb potential',
                                                 'ngridq':'q-point grid sizes','ngridk':'k-point grid size',
                                                 'broad':'Lorentzian broadening for all spectra (in hartree)',
                                                 'screentype':"""Defines which type of screening is to be used.
Type:   choose from:
full
diag
noinvdiag
longrange
""",
                                                 'nempty':"""Number of empty states. This parameter determines the energy cutoff for the excitation spectra. For determining the number of states related to an energy cutoff, perform one iteration of a SCF calculation, setting nempty to a higher value and check the EIGVAL.OUT.
Type: 	integer
"""  ,
                                                 'xstype':"""Attribute: xstype

Should TDDFT be used or BSE.
Type: 	choose from:
TDDFT
BSE""",
                                                 'use gw':'use gw band structure as basis (true or false)',
                                                 'bsetype':"""Defines which parts of the BSE Hamiltonian are to be considered.
Type: 	choose from:
IP (Independent particles)
RPA (random phase approximation)
singlet
triplet"""}

        self.relax_file_timestamp = None

    def find_engine_folder(self):
        p = subprocess.Popen([search_command, 'excitingser'], stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=shell_bool)
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

    def start_ground_state(self, crystal_structure, band_structure_points=None,blocking=False):
        """This method starts a ground state calculation in a subprocess. The configuration is stored in scf_options.

Args:
    - crystal_structure:        A CrystalStructure or MolecularStructure object that represents the geometric structure of the material under study.

Keyword Args:
    - band_structure_points:    If supplied automatically triggers a calculation of the band structure of the material.
                                Must be a list of two component lists as following:
                                [[np.array([0,0,0]),'Gamma'],[np.array([0.5,0.5,0.5]),'W']]
                                Default: None

    - blocking:                 Determines whether the function call will block the main process or run in the background.
                                Helpful when looping over different calculations in the builtin python terminal.
                                Default: False
Returns:
    - None
        """
        try:
            os.remove(self.project_directory + self.working_dirctory + '/INFO.OUT')
        except Exception as e:
            print(e)
        self._read_timestamps()
        self.current_output_file = 'INFO.OUT'


        tree = self._make_tree()
        self._add_scf_to_tree(tree, crystal_structure)
        if band_structure_points is not None:
            self._add_bs_to_tree(tree, band_structure_points)
        self._write_input_file(tree)
        time.sleep(0.05)
        self._start_engine(blocking=blocking)

    def start_optical_spectrum(self,crystal_structure):
        """This method starts a optical spectrum calculation in a subprocess. The configuration is stored in optical_spectrum_options.

Args:
    - crystal_structure:        A CrystalStructure or MolecularStructure object that represents the geometric structure of the material under study.

Returns:
    - None
        """

        self._filenames_tasks['optical spectrum'] = '/EPSILON_BSE' + self.optical_spectrum_options['bsetype'] + '_SCRfull_OC11.OUT'
        self._read_timestamps()
        self.current_output_file = 'INFOXS.OUT'

        tree = self._make_tree()
        self._add_scf_to_tree(tree, crystal_structure,skip=True)
        self._add_optical_spectrum_to_tree(tree)
        self._write_input_file(tree)
        time.sleep(0.05)
        self._start_engine()

    def start_gw(self,crystal_structure,band_structure_points=None,blocking=False):
        """This method starts a g0w0 calculation in a subprocess. The configuration is stored in gw_options.

Args:
    - crystal_structure:        A CrystalStructure or MolecularStructure object that represents the geometric structure of the material under study.

Keyword Args:
    - band_structure_points:    If supplied automatically triggers a calculation of the band structure of the material.
                                Must be a list of two component lists as following:
                                [[np.array([0,0,0]),'Gamma'],[np.array([0.5,0.5,0.5]),'W']]
                                Default: None

    - blocking:                 Determines whether the function call will block the main process or run in the background.
                                Helpful when looping over different calculations in the builtin python terminal.
                                Default: False
Returns:
    - None
        """
        self._read_timestamps()
        self.current_output_file = 'GW_INFO.OUT'
        tree = self._make_tree()
        self._add_scf_to_tree(tree, crystal_structure,skip=True)
        self._add_gw_to_tree(tree, taskname='g0w0')
        self._write_input_file(tree)

        def first_round():
            self._start_engine(blocking=blocking)
            if self.custom_command_active:
                while self.is_engine_running(tasks=['g0w0']):
                    time.sleep(1)
            else:
                self.engine_process.wait()
            if band_structure_points is not None:
                tree = self._make_tree()
                self._add_scf_to_tree(tree, crystal_structure,skip=True)
                self._add_gw_to_tree(tree, taskname='band')
                self._add_bs_to_tree(tree, band_structure_points)
                self._write_input_file(tree)
                self._start_engine(blocking=blocking)

        t = threading.Thread(target=first_round)
        t.start()

    def start_phonon(self, crystal_structure, band_structure_points):
        """This method starts a phonon bandstructure calculation in a subprocess. The configuration is stored in phonons_options.

Args:
    - crystal_structure:        A CrystalStructure or MolecularStructure object that represents the geometric structure of the material under study.

    - band_structure_points:    If supplied automatically triggers a calculation of the band structure of the material.
                                Must be a list of two component lists as following:
                                [[np.array([0,0,0]),'Gamma'],[np.array([0.5,0.5,0.5]),'W']]
                                Default: None

Returns:
    - None
        """
        self._read_timestamps()
        tree = self._make_tree()
        self._add_scf_to_tree(tree, crystal_structure)
        self._add_phonon_to_tree(tree, band_structure_points)
        self._write_input_file(tree)
        time.sleep(0.05)
        self._start_engine()

    def start_relax(self,crystal_structure):
        """This method starts a structure relaxation calculation in a subprocess. The configuration is stored in relax_options.

Args:
    - crystal_structure:        A CrystalStructure or MolecularStructure object that represents the geometric structure of the material under study.

Returns:
    - None
        """
        self._read_timestamps()
        tree = self._make_tree()
        self._add_scf_to_tree(tree, crystal_structure)
        self._add_relax_to_tree(tree)
        self._write_input_file(tree)
        time.sleep(0.05)
        self._start_engine()

    def load_relax_structure(self):
        """This method loads the result of a relaxation calculation, which is a molecular or crystal structure.

Returns:
    - CrystalStructure or MolecularStructure object depending on the material under study.
        """
        file = self.project_directory+self.working_dirctory + 'geometry_opt.xml'
        if not os.path.isfile(file):
            return None
        if self.relax_file_timestamp is not None and os.path.getmtime(file) == self.relax_file_timestamp:
            return None
        struc = self.parse_input_file(file)
        self.relax_file_timestamp = os.path.getmtime(file)
        return struc

    def read_scf_status(self):
        """This method reads the result of a self consistent ground state calculation.

Returns:
    - res: Nx2 numpy array with iteration number and scf energy in the first and second column respectively.
        """
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

    def read_bandstructure(self,crystal_structure=None,special_k_points=None):
        """This method reads the result of a electronic band structure calculation.

Keyword args:
    - crystal_structure:    For some engines the crystal structure must be re-supplied.
                            Default: None

    - special_k_points:     For some engines the special k-points must be re-supplied
                            Default: None

Returns:
    - band_structure:       A BandStructure object with the latest band structure result found.
        """
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

        special_k_points_together = zip(special_k_points, special_k_points_label)


        return sst.BandStructure(bands, special_k_points=special_k_points_together)

    def read_gw_bandstructure(self,filename='BAND-QP.OUT',special_k_points=None,structure=None):
        """This method reads the result of a gw electronic band structure calculation.

Keyword args:
    - filename:             filename to be read.

Returns:
    - band_structure:       A BandStructure object with the latest band structure result found.
        """

        if special_k_points is None or structure is None:
            band_structure = self.read_bandstructure()
            special_k_points = band_structure.special_k_points
        else:
            special_k_points = sst.calculate_path_length(structure,special_k_points)



        band_qp = np.loadtxt(self.project_directory + self.working_dirctory + filename)
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

        if filename == 'PHDISP.OUT':
            bs_type = 'phonon'
        else:
            bs_type = 'electronic'
        return sst.BandStructure(bands_qp,special_k_points=special_k_points,bs_type=bs_type)

    def read_phonon_bandstructure(self,special_k_points=None,structure=None):
        """This method reads the result of a phonon band structure calculation.

Returns:
    - band_structure:       A BandStructure object with the latest phonon band structure result found.
        """
        return self.read_gw_bandstructure(filename='PHDISP.OUT',special_k_points=special_k_points,structure=structure)

    def read_optical_spectrum(self):
        """This method reads the result of a optical spectrum calculation.

Returns:
    - optical_spectrum:       A OpticalSpectrum object with the latest optical spectrum result found.
        """
        eps_11 = np.loadtxt(self.project_directory + self.working_dirctory + 'EPSILON_BSE' + self.optical_spectrum_options['bsetype'] + '_SCRfull_OC11.OUT')
        eps_22 = np.loadtxt(self.project_directory + self.working_dirctory + 'EPSILON_BSE' + self.optical_spectrum_options['bsetype'] + '_SCRfull_OC22.OUT')
        eps_33 = np.loadtxt(self.project_directory + self.working_dirctory + 'EPSILON_BSE' + self.optical_spectrum_options['bsetype'] + '_SCRfull_OC33.OUT')


        list_of_eps2 = [eps_11[:,2],eps_22[:,2],eps_33[:,2]]
        list_of_eps1 = [eps_11[:,1],eps_22[:,1],eps_33[:,1]]

        return sst.OpticalSpectrum(eps_11[:,0]*hartree,list_of_eps2,epsilon1=list_of_eps1)

    def read_ks_state(self):
        """This method reads the result of a electronic state calculation and returns the modulo squared of the wavefunction.

Returns:
    - ks_density:       A KohnShamDensity or MolecularDensity object with the latest result found.
        """
        self._convert_3d_plot()
        l_data = np.genfromtxt(self.project_directory + self.working_dirctory + 'WF3D.xsf', skip_header=9, skip_footer=2, dtype=np.float)
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
        """This method starts a calculation of a specific electronic state in a subprocess.

Args:
    - crystal_structure:        A CrystalStructure or MolecularStructure object that represents the geometric structure of the material under study.

    - bs_point:                 Index of the k point and band of the state that should be calculated. Should be a list with [k,n].
                                In case of a molcular calculation the k index must still be supplied but is not used.

Keyword Args:
    - grid:                     Determines the spatial grid to be used. Must be a string that is valid for the respective engine.
                                Default: '40 40 40'

Returns:
    - None
                """
        self._read_timestamps()
        tree = self._make_tree()
        self._add_scf_to_tree(tree, crystal_structure, skip=True)
        self._add_ks_density_to_tree(tree, bs_point, grid)
        self._write_input_file(tree)
        self._start_engine()

    def calculate_electron_density(self, crystal_structure, grid='40 40 40'):
        """This method starts a calculation of the total (pseudo-) electron density in a subprocess.

Args:
    - crystal_structure:        A CrystalStructure or MolecularStructure object that represents the geometric structure of the material under study.

Keyword Args:
    - grid:                     Determines the spatial grid to be used. Must be a string that is valid for the respective engine.
                                Default: '40 40 40'

Returns:
    - None
                """
        # TODO the total energy density calculation for exciting (if possible)
        raise NotImplementedError

    def kill_engine(self):
        """Stops the execution of the engine process. Only possible for local execution and not in case of cluster calculation"""
        try:
            self.engine_process.kill()
        except Exception as e:
            print(e)

    def is_engine_running(self,tasks=None):
        """Determines whether the engine is currently running.

Keyword args:
    - tasks:    List of tasks that are supposed to be running. Must be supplied when the calculations run on a cluster.
                Possible tasks are: ['bandstructure', 'relax', 'ks density', 'scf', 'g0w0', 'g0w0 bands', 'optical spectrum']

Returns:
    - res:      Boolean result. True: engine is running. False: engine is not running.
"""
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

    def reset_to_defaults(self):
        """Reset all configurations to their defaults."""
        default_handler = Handler()
        self.scf_options.update(default_handler.scf_options)
        self.gw_options.update(default_handler.gw_options)
        self.optical_spectrum_options.update(default_handler.optical_spectrum_options)
        self.relax_options.update(default_handler.relax_options)
        self.phonons_options.update(default_handler.phonons_options)

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
                                  nempty=self.optical_spectrum_options['nempty_screening'])
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
            if filename and os.path.isfile(self.project_directory+self.working_dirctory+filename):
                self._timestamp_tasks[task]=os.path.getmtime(self.project_directory + self.working_dirctory + filename)

    def _check_if_scf_is_finished(self):
        try:
            f = open(self.project_directory + self.working_dirctory + self.info_file, 'r')
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
        if not os.path.isdir(self.project_directory + self.working_dirctory):
            os.mkdir(self.project_directory + self.working_dirctory)
        xmlstr = minidom.parseString(ET.tostring(tree.getroot())).toprettyxml(indent="   ")
        with open(self.project_directory + self.working_dirctory + self.input_filename, "w") as f:
            f.write(xmlstr)

            # tree.write(self.project_directory+self.working_dirctory+self.input_filename)

    def _start_engine(self,blocking=False):
        os.chdir(self.project_directory + self.working_dirctory)
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
        if blocking:
            while self.is_engine_running():
                time.sleep(0.1)
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
            file_exists = os.path.isfile(self.project_directory + self.working_dirctory + filename)
            # files_exists[task] = file_exists
            if file_exists:
                timestamp = os.path.getmtime(self.project_directory + self.working_dirctory + filename)
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
        os.chdir(self.project_directory + self.working_dirctory)
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
    handler.project_directory = "/home/jannick/OpenDFT_projects/test_abinit"

    trash_bs_points = np.array([[0, 0, 0], [0.50, 0.500, 0.0], [0.500, 0.500, 0.500]
                                   , [0.000, 0.000, 0.000], [0.500, 0.500, 0.000], [0.750, 0.500, 0.250],
                                [0.750, 0.375, 0.375], [0.000, 0.000, 0.000]])
    trash_bs_labels = ['GAMMA', 'X', 'L', 'GAMMA', 'X', 'W', 'K', 'GAMMA']
    path = list(zip(trash_bs_points, trash_bs_labels))

    atoms = np.array([[0, 0, 0, 6], [0.25, 0.25, 0.25, 6]])
    unit_cell = 6.6 * np.array([[0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])

    crystal_structure = sst.CrystalStructure(unit_cell, atoms)
    handler.calculate_ks_density(crystal_structure,[1,1])

    handler.read_gw_bandstructure(special_k_points=path,structure=crystal_structure)


    # handler.custom_command_active = True
    # ac_bo  = handler.is_engine_running(tasks = ['scf','bandstructure','g0w0'])
    # ac_bo2 = handler._check_if_scf_is_finished()

    # print(ac_bo,ac_bo2)

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


