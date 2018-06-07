import os
import subprocess
from abinit_handler import Handler as AbinitHandler
from abinit_handler import p_table,p_table_rev
from quantum_espresso_handler import Handler as QeHandler
from little_helpers import convert_to_ordered,find_data_file
import solid_state_tools as sst
import numpy as np
import shutil
from sortedcontainers import SortedSet
import scipy.special



class OceanHandler(object):

    def __init__(self):
        self.supported_methods.add('optical spectrum')
        self.working_dirctory = '/abinit_ocean_files/'

        self.info_text = self.info_text + """
        
ocean is  an ab initio Density  Functional  Theory  (DFT)  +  Bethe-Salpeter  Equation  (BSE)
code for calculations of core-level spectra.  Currently the code allows for the calculations of x-
ray absorption spectra (XAS), x-ray emission (XES), non-resonant x-ray inelastic x-ray spectra
(NRIXS), and (direct) resonant inelastic x-ray scattering (RIXS) of periodic systems.  The code
is written in Fortran 90 with associated shell and Perl scripting."""

        self.optical_spectrum_options = convert_to_ordered({'diemac':'5.0','CNBSE.xmesh':'6 6 6','screen.shells':'4.0','screen.shells':'3.5','cnbse.rad':'3.5',
        'cnbse.broaden':'0.1','edges':'0 1 0','screen.nbands':'80','screen.nkpt':'2 2 2','opts.core_states':'1 0 0 0',
                                                            'opts.relativity':'scalar rel','opts.functional':'lda','opts.valence':'2.0 3.5 0.0 0.0',
        'fill.pow':'2','fill.energy':'0.30 2.00 0.0001','fill.cutoff':'3.5','fill.fourier':'0.05 20',
        'operator':'dipole','polarization':'all','k vector':'0 1 0','photon energy':'600'})

    def start_optical_spectrum(self, crystal_structure):
        """This method starts a optical spectrum calculation in a subprocess. The configuration is stored in optical_spectrum_options.

Args:
    - crystal_structure:        A CrystalStructure or MolecularStructure object that represents the geometric structure of the material under study.

Returns:
    - None
        """

        self.current_input_file = 'ocean.in'
        self.current_output_file = 'ocean.out'

        f = self._make_ocean_input_file()
        self._add_ocean_to_file(f,crystal_structure)
        f.close()

        self._make_opts_file(crystal_structure)
        self._make_fill_file()
        self._make_photon_file()

        self._copy_default_ocean_pseudos(crystal_structure)

        if os.path.isdir(self.project_directory + self.working_dirctory+'/CNBSE'):
            shutil.rmtree(self.project_directory + self.working_dirctory+'/CNBSE')

        self._start_ocean_engine()

    def read_optical_spectrum(self):
        """This method reads the result of a optical spectrum calculation.

Returns:
    - optical_spectrum:       A OpticalSpectrum object with the latest optical spectrum result found.
        """

        def load_eps(i):
            files = []
            for file in os.listdir(self.project_directory+self.working_dirctory+'/CNBSE'):
                if file.startswith("absspct_") and file.endswith('_0{}'.format(i+1)):
                    files.append(file)

            if not files:
                return None,None,None


            data = np.loadtxt(self.project_directory+self.working_dirctory+'/CNBSE/'+files[0], skiprows=2)

            for file in files[1:]:
                data += np.loadtxt(self.project_directory+self.working_dirctory+'/CNBSE/'+file, skiprows=2)

            energy = data[:, 0]
            eps2 = data[:, 1]
            eps1 = data[:, 3]
            return energy,eps1,eps2

        eps1_list = []
        eps2_list = []
        for i in range(3):
            energy,eps1,eps2 = load_eps(i)
            eps1_list.append(eps1)
            eps2_list.append(eps2)

        return sst.OpticalSpectrum(energy, eps2_list, epsilon1=eps1_list)

    def _make_ocean_input_file(self, filename='ocean.in'):
        if not os.path.isdir(self.project_directory + self.working_dirctory):
            os.mkdir(self.project_directory + self.working_dirctory)
        f = open(self.project_directory + self.working_dirctory + '/' + filename, 'w')
        return f

    def _add_ocean_to_file(self,file,crystal_structure):
        file.write('dft{ abi }\n')
        file.write(r"ppdir {'../'}"+'\n')

        self._add_scf_to_file(file,crystal_structure,brakets=True)

        file.write('nbands '+self.scf_options['nband']+'\n')

        atom_species = SortedSet(crystal_structure.atoms[:,3])
        file.write('pp_list{\n')
        for specie in atom_species:
            file.write(p_table[specie]+'.fhi\n' )
        file.write('}\n')

        file.write('\n# BSE options \n')
        ocean_opts = dict((k, self.optical_spectrum_options[k]) for k in
                          ('diemac','CNBSE.xmesh','screen.shells','cnbse.rad','cnbse.broaden','edges','screen.nbands','screen.nkpt') if k in self.optical_spectrum_options)
        for key, item in ocean_opts.items():
            file.write(key+'{ '+item+' }\n')

        edge_0 = int(self.optical_spectrum_options['edges'].split()[0])
        if edge_0 < 0:
            Z = abs(edge_0)
        else:
            Z = SortedSet(crystal_structure.atoms[:,3].astype('int'))[edge_0-1]

        file.write("""
opf.opts{{ {0} ocean.opts }}
opf.fill{{ {0} ocean.fill }}""".format(Z))

    def _start_ocean_engine(self,blocking=False):
        os.chdir(self.project_directory + self.working_dirctory)
        if self.custom_command_active:
            command = ['bash', self.custom_command]
        else:
            command = self._engine_command


        self.engine_process = subprocess.Popen('exec ocean.pl ocean.in >ocean.out', stdout=subprocess.PIPE,
                                               stderr=subprocess.PIPE,shell=True)
        os.chdir(self.project_directory)
        if blocking:
            while self.is_engine_running():
                time.sleep(0.1)

    def _make_opts_file(self,crystal_structure,filename='ocean.opts'):
        if not os.path.isdir(self.project_directory + self.working_dirctory):
            os.mkdir(self.project_directory + self.working_dirctory)
        f = open(self.project_directory + self.working_dirctory + '/' + filename, 'w')

        edge_0 = int(self.optical_spectrum_options['edges'].split()[0])
        if edge_0 < 0:
            Z = abs(edge_0)
        else:
            Z = SortedSet(crystal_structure.atoms[:,3].astype('int'))[edge_0-1]

        f.write('{0:03d}\n'.format(Z))
        f.write(self.optical_spectrum_options['opts.core_states'] +'\n')
        f.write(self.optical_spectrum_options['opts.relativity'] + '\n')
        f.write(self.optical_spectrum_options['opts.functional'] + '\n')
        f.write(self.optical_spectrum_options['opts.valence'] + '\n')
        f.write(self.optical_spectrum_options['opts.valence'] + '\n')

        f.close()

    def _make_fill_file(self,filename='ocean.fill'):
        if not os.path.isdir(self.project_directory + self.working_dirctory):
            os.mkdir(self.project_directory + self.working_dirctory)
        f = open(self.project_directory + self.working_dirctory + '/' + filename, 'w')
        f.write(self.optical_spectrum_options['fill.pow']+'\n')
        f.write(self.optical_spectrum_options['fill.energy']+'\n')
        f.write(self.optical_spectrum_options['fill.cutoff']+'\n')
        f.write(self.optical_spectrum_options['fill.fourier']+'\n')

        f.close()

    def _make_photon_file(self):
        if not os.path.isdir(self.project_directory + self.working_dirctory):
            os.mkdir(self.project_directory + self.working_dirctory)

        if self.optical_spectrum_options['polarization'].strip().lower() == 'all':
            polarizations = ['1 0 0','0 1 0','0 0 1']
        else:
            polarizations = [self.optical_spectrum_options['polarization']]

        for i,polarization in enumerate(polarizations):
            f = open(self.project_directory + self.working_dirctory + '/' + 'photon{0}'.format(i+1), 'w')

            f.write(self.optical_spectrum_options['operator']+'\n')
            f.write('cartesian '+polarization+'\n')
            f.write('end\n')
            f.write('cartesian '+self.optical_spectrum_options['k vector']+'\n')
            f.write('end\n')
            f.write(self.optical_spectrum_options['photon energy'])

        f.close()

    def _copy_default_ocean_pseudos(self, crystal_structure):
        atoms = sorted(set(crystal_structure.atoms[:,3]))
        atoms_names = [p_table[atom] for atom in atoms]
        installation_folder = find_data_file('')

        if not os.path.isdir(self.project_directory+self.pseudo_directory):
            os.mkdir(self.project_directory+self.pseudo_directory)

        pseudo_files = []
        for atom in atoms_names:
            file = atom+'.fhi'
            pseudo_files.append(file)
            filepath = self.project_directory+self.working_dirctory+file
            if not os.path.isfile(filepath):
                shutil.copyfile(installation_folder+'/data/pseudos/ocean/'+file,filepath)

        return pseudo_files


class FeffHandler(object):
    def __init__(self):
        self.supported_methods.add('optical spectrum')
        self.working_dirctory = '/abinit_feff_files/'

        self.info_text = self.info_text + """ """

        self.optical_spectrum_options = convert_to_ordered(
            {'sphere radius':'6.0','atom':'0','temperature':'300','debye temperature':'500','edge energy':'9000','edge amplitude':'1.0','edge width':'10'})

    def start_optical_spectrum(self, crystal_structure):
        """This method starts a optical spectrum calculation in a subprocess. The configuration is stored in optical_spectrum_options.

Args:
    - crystal_structure:        A CrystalStructure or MolecularStructure object that represents the geometric structure of the material under study.

Returns:
    - None
        """

        self.current_input_file = 'feff.inp'
        self.current_output_file = 'feff.out'

        f = self._make_feff_input_file()
        self._add_feff_to_file(f,crystal_structure)
        f.close()

        self._start_feff_engine()

    def _make_feff_input_file(self, filename='feff.inp'):
        if not os.path.isdir(self.project_directory + self.working_dirctory):
            os.mkdir(self.project_directory + self.working_dirctory)
        f = open(self.project_directory + self.working_dirctory + '/' + filename, 'w')
        return f

    def _add_feff_to_file(self,file,crystal_structure):
        file.write('DEBYE {0} {1}  \n\n'.format(self.optical_spectrum_options['temperature'],self.optical_spectrum_options['debye temperature']))

        scattering_atom = int(self.optical_spectrum_options['atom'])
        sphere_radius = float(self.optical_spectrum_options['sphere radius'])

        single_cell_coord = crystal_structure.calc_absolute_coordinates()
        Z_scattering = int( single_cell_coord[scattering_atom,3] )

        atoms = self._find_atoms_within_sphere(crystal_structure,sphere_radius,scattering_atom)

        species = SortedSet(atoms[:,3].astype('int'))

        file.write('POTENTIALS\n')
        file.write('  0 {}\n'.format(Z_scattering))
        for i,specie in enumerate(species):
            file.write('  {0} {1}\n'.format(i+1,specie))

        file.write('\nATOMS\n')
        for atom in atoms:
            coords = atom[:3]
            Z = int(atom[3])
            if np.linalg.norm(coords-single_cell_coord[scattering_atom,:3]) <1e-6:
                potential_number = 0
            else:
                potential_number = species.index(Z)+1

            in_list = [coords[0]*sst.bohr,coords[1]*sst.bohr,coords[2]*sst.bohr,potential_number]
            file.write('  {0:1.10f} {1:1.10f} {2:1.10f} {3}\n'.format(*in_list))

    def _start_feff_engine(self,blocking=False):
        os.chdir(self.project_directory + self.working_dirctory)
        if self.custom_command_active:
            command = ['bash', self.custom_command]
        else:
            command = self._engine_command


        self.engine_process = subprocess.Popen('exec feff.x >feff.out', stdout=subprocess.PIPE,
                                               stderr=subprocess.PIPE,shell=True)
        os.chdir(self.project_directory)
        if blocking:
            while self.is_engine_running():
                time.sleep(0.1)

    def _find_atoms_within_sphere(self,crystal_structure,sphere_radius,atom):
        atomic_coord_single = crystal_structure.calc_absolute_coordinates()
        atomic_coord_rep = crystal_structure.calc_absolute_coordinates(repeat=(5,5,5),offset=(2,2,2))
        dist = np.linalg.norm(atomic_coord_rep[:,:3] - atomic_coord_single[atom,:3],axis=1)
        sphere_mask = dist < sphere_radius
        sphere_atoms = atomic_coord_rep[sphere_mask,:]
        return sphere_atoms

    def read_optical_spectrum(self):

        try:
            f = open(self.project_directory + self.working_dirctory + u'/chi.dat', 'r')
        except IOError:
            return None
        text = f.read()
        f.close()

        lines = text.split('\n')
        for i,line in enumerate(lines):
            if line.strip().startswith('0.00'):
                break

        data = np.loadtxt(self.project_directory + self.working_dirctory + u'/chi.dat',skiprows=i)

        Ek = float(self.optical_spectrum_options['edge energy'])
        hbar = 1.0545718e-34
        me = 9.10938356e-31
        a0 = 5.2917721092e-11
        e = 1.60217662e-19
        E = (data[:, 0] / 1e-10) ** 2 * hbar ** 2 / (2 * me) / e + Ek

        width = float(self.optical_spectrum_options['edge width'])
        A0 = float(self.optical_spectrum_options['edge amplitude'])
        edge_f = lambda E: (scipy.special.erf((E - Ek) / width) + 1)/2 / (E / Ek) ** 4

        region = E.max()-E.min()

        E_plot = np.linspace(E.min()-region*0.1,E.max(),2000)


        edge = edge_f(E_plot)
        exaf_interp = np.interp(E_plot, E, data[:, 1], left=0)

        res_abs = A0*(exaf_interp + edge)
        # c = 299792458
        # omega = E_plot*e/hbar
        # eps2 = res_abs*1e6*c/omega

        return sst.OpticalSpectrum(E_plot,res_abs)


class OceanAbinit(OceanHandler,AbinitHandler):

    def __init__(self):
        AbinitHandler.__init__(self)
        OceanHandler.__init__(self)


class OceanQe(OceanHandler,QeHandler):

    def __init__(self):
        QeHandler.__init__(self)
        OceanHandler.__init__(self)


class FeffAbinit(FeffHandler,AbinitHandler):

    def __init__(self):
        AbinitHandler.__init__(self)
        FeffHandler.__init__(self)


if __name__ == '__main__':
    # ocean_abi_handler = OceanAbinit()
    # ocean_qe_handler = OceanQe()
    #
    # ocean_abi_handler.project_directory = '/home/jannick/OpenDFT_projects/ocean_test'
    #
    # ocean_abi_handler.optical_spectrum_options['edges'] = '-6 1 0'
    # ocean_abi_handler.optical_spectrum_options['opts.valence'] = '0.8 2.0 0.2 0.0'

    atoms = np.array([[0, 0, 0, 6], [0.25, 0.25, 0.25, 6]])
    unit_cell = 6.719 * np.array([[0.0,0.5,0.5], [0.5, 0, 0.5], [0.5, 0.5, 0.0]])

    crystal_structure = sst.CrystalStructure(unit_cell, atoms)

    # ocean_abi_handler.start_optical_spectrum(crystal_structure)

    # spec = ocean_abi_handler.read_optical_spectrum()

    # import matplotlib.pyplot as plt
    #
    # plt.plot(spec.energy,spec.epsilon2)
    handler = FeffAbinit()
    handler.project_directory = '/home/jannick/OpenDFT_projects/feff_testing'
    handler.optical_spectrum_options['sphere radius'] = '9.0'
    handler.start_optical_spectrum(crystal_structure)
    import time
    time.sleep(1.0)
    spec = handler.read_optical_spectrum()

    import matplotlib.pyplot as plt

    plt.plot(spec.energy,spec.epsilon2)