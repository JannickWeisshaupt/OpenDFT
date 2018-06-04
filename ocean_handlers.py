import os
import subprocess
from abinit_handler import Handler as AbinitHandler
from abinit_handler import p_table,p_table_rev
from quantum_espresso_handler import Handler as QeHandler
from little_helpers import convert_to_ordered
import solid_state_tools as sst
import numpy as np
import shutil
from sortedcontainers import SortedSet


class OceanHandler(object):

    def __init__(self):
        self.supported_methods.add('optical spectrum')
        self.working_dirctory = '/abinit_ocean_files/'

        self.optical_spectrum_options = convert_to_ordered({'diemac':'5.0','CNBSE.xmesh':'6 6 6','screen.shells':'4.0','screen.shells':'3.5','cnbse.rad':'3.5',
        'cnbse.broaden':'0.1','edges':'0 1 0','screen.nbands':'80','screen.nkpt':'2 2 2','opts.core_states':'1 0 0 0',
                                                            'opts.relativity':'scalar rel','opts.functional':'lda','opts.valence':'2.0 3.5 0.0 0.0',
        'fill.pow':'2','fill.energy':'0.30 2.00 0.0001','fill.cutoff':'3.5','fill.fourier':'0.05 20',
        'operator':'dipole','polarization':'all','k vector':'0 1 0','energy range':'600'})

    def start_optical_spectrum(self, crystal_structure):
        """This method starts a optical spectrum calculation in a subprocess. The configuration is stored in optical_spectrum_options.

Args:
    - crystal_structure:        A CrystalStructure or MolecularStructure object that represents the geometric structure of the material under study.

Returns:
    - None
        """

        if os.path.isdir(self.project_directory + self.working_dirctory+'/CNBSE'):
            shutil.rmtree(self.project_directory + self.working_dirctory+'/CNBSE')

        self.current_input_file = 'ocean.in'
        self.current_output_file = 'ocean.out'

        f = self._make_ocean_input_file()
        self._add_ocean_to_file(f,crystal_structure)
        f.close()

        self._make_opts_file()
        self._make_fill_file()
        self._make_photon_file()

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

        file.write("""
opf.opts{{ {0} ocean.opts }}
opf.fill{{ {0} ocean.fill }}""".format(abs(int(self.optical_spectrum_options['edges'].split()[0]) )))

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

    def _make_opts_file(self,filename='ocean.opts'):
        if not os.path.isdir(self.project_directory + self.working_dirctory):
            os.mkdir(self.project_directory + self.working_dirctory)
        f = open(self.project_directory + self.working_dirctory + '/' + filename, 'w')
        Z = abs(int(self.optical_spectrum_options['edges'].split()[0]))
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
            f.write(self.optical_spectrum_options['energy range'])

        f.close()

class OceanAbinit(OceanHandler,AbinitHandler):

    def __init__(self):
        AbinitHandler.__init__(self)
        OceanHandler.__init__(self)


class OceanQe(OceanHandler,QeHandler):

    def __init__(self):
        QeHandler.__init__(self)
        OceanHandler.__init__(self)


if __name__ == '__main__':
    ocean_abi_handler = OceanAbinit()
    ocean_qe_handler = OceanQe()

    ocean_abi_handler.project_directory = '/home/jannick/OpenDFT_projects/ocean_test'

    ocean_abi_handler.optical_spectrum_options['edges'] = '-6 1 0'
    ocean_abi_handler.optical_spectrum_options['opts.valence'] = '0.8 2.0 0.2 0.0'

    atoms = np.array([[0, 0, 0, 6], [0.25, 0.25, 0.25, 6]])
    unit_cell = 6.719 * np.array([[0.0,0.5,0.5], [0.5, 0, 0.5], [0.5, 0.5, 0.0]])

    crystal_structure = sst.CrystalStructure(unit_cell, atoms)

    ocean_abi_handler.start_optical_spectrum(crystal_structure)

    # spec = ocean_abi_handler.read_optical_spectrum()

    # import matplotlib.pyplot as plt
    #
    # plt.plot(spec.energy,spec.epsilon2)
