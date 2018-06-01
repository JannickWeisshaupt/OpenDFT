import os
import subprocess
from abinit_handler import Handler as AbinitHandler
from quantum_espresso_handler import Handler as QeHandler
from little_helpers import convert_to_ordered
import solid_state_tools as sst
import numpy as np




class OceanHandler(object):

    def __init__(self):
        self.supported_methods.add('optical spectrum')

        self.optical_spectrum_options = convert_to_ordered({'diemac':'5.0','CNBSE.xmesh':'6 6 6','screen.shells':'4.0'})


    def start_optical_spectrum(self, crystal_structure):
        """This method starts a optical spectrum calculation in a subprocess. The configuration is stored in optical_spectrum_options.

Args:
    - crystal_structure:        A CrystalStructure or MolecularStructure object that represents the geometric structure of the material under study.

Returns:
    - None
        """

        # f = self._make_ocean_input_file()
        # self._add_ocean_to_file(f,crystal_structure)
        # f.close()
        self._start_ocean_engine()

    def _make_ocean_input_file(self, filename='ocean.in'):
        if not os.path.isdir(self.project_directory + self.working_dirctory):
            os.mkdir(self.project_directory + self.working_dirctory)
        f = open(self.project_directory + self.working_dirctory + '/' + filename, 'w')
        return f

    def _add_ocean_to_file(self,file,crystal_structure):
        pass

    def _start_ocean_engine(self,blocking=False):
        os.chdir(self.project_directory + self.working_dirctory)
        if self.custom_command_active:
            command = ['bash', self.custom_command]
        else:
            command = self._engine_command


        self.engine_process = subprocess.Popen(['ocean.pl','ocean.in'], stdout=subprocess.PIPE,
                                               stderr=subprocess.PIPE)
        os.chdir(self.project_directory)
        if blocking:
            while self.is_engine_running():
                time.sleep(0.1)



    def read_optical_spectrum(self):
        """This method reads the result of a optical spectrum calculation.

Returns:
    - optical_spectrum:       A OpticalSpectrum object with the latest optical spectrum result found.
        """
        data = np.loadtxt(
            self.project_directory + self.working_dirctory + '/CNBSE/absspct_C_.0001_1s_01',skiprows=2)

        energy = data[:,0]
        eps2 = data[:,1]
        eps1 = data[:,3]


        return sst.OpticalSpectrum(energy,eps2,epsilon1=eps1)

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
    # ocean_abi_handler.start_optical_spectrum(None)
    spec = ocean_abi_handler.read_optical_spectrum()

    import matplotlib.pyplot as plt

    plt.plot(spec.energy,spec.epsilon2)
