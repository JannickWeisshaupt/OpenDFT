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
from six import string_types

atomic_mass = pt.mass
p_table = {i: el.__repr__() for i, el in enumerate(pt.elements)}
p_table_rev = {el.__repr__(): i for i, el in enumerate(pt.elements)}
atomic_mass = {i: el.mass for i, el in enumerate(pt.elements)}

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
        self._engine_command = ["pw.x"]
        self._working_dirctory = '/quantum_espresso_files/'
        self.pseudo_directory = '/pseudos/'
        self.engine_process = None
        self._info_file = 'scf.out'
        self._filenames_tasks = {}

        self._timestamp_tasks = {}

        self.project_directory = None
        self._input_filename = 'input.xml'
        self.custom_command = ''
        self.custom_command_active = False
        self.exciting_folder = self.find_engine_folder()
        self.scf_options = {'prefix':'title','ecutwfc':'12.0','k points':'6 6 6','k point shift':'1 1 1','k points band':'30','nbnd':'10'}
        # self.scf_options_non_string_type  = {'ecutwfc':float}
        self.scf_options_tooltip = {}

        self.general_options = {'title': 'title'}
        self.bs_options = {}
        self.relax_options = {}
        self.gw_options = {}
        self.gw_options_tooltip = {}
        self.phonons_options = {}
        self.phonons_options_tooltip = {}
        self.optical_spectrum_options = {}

    def find_engine_folder(self):
        p = subprocess.Popen(['which', 'pw.x'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        res, err = p.communicate()
        res = res.decode()
        res = res.split('bin')[0]
        return res.strip()

    def parse_input_file(self, filename):
        raise NotImplementedError

    def start_ground_state(self, crystal_structure, band_structure_points=None):
        # self._correct_types()
        file = self._make_input_file()
        self._add_scf_to_file(file,crystal_structure)
        file.close()
        self._start_engine()

        if band_structure_points is not None:
            def run_bs():
                while self.is_engine_running():
                    time.sleep(0.5)
                file = self._make_input_file(filename='bands.in')
                self._add_scf_to_file(file,crystal_structure,calculation='bands',band_points=band_structure_points)
                file.close()
                self._start_engine(filename='bands.in')


            t = threading.Thread(target=run_bs)
            t.start()

    def start_optical_spectrum(self, crystal_structure):
        raise NotImplementedError

    def start_gw(self, crystal_structure, band_structure_points=None):
        raise NotImplementedError

    def start_phonon(self, crystal_structure, band_structure_points):
        raise NotImplementedError

    def start_relax(self, crystal_structure):
        raise NotImplementedError

    def load_relax_structure(self):
        raise NotImplementedError

    def read_scf_status(self):
        try:
            f = open(self.project_directory + self._working_dirctory + self._info_file, 'r')
        except IOError:
            return None
        info_text = f.read()
        f.close()


        scf_energy_list = []
        matches = re.findall(r"total energy[\s\t]*=[\s\t]*[-+]?\d*\.\d+", info_text)
        for match in matches:
            ms = match.split('=')
            scf_energy_list.append(float(ms[1]))

        res = np.array(zip(range(len(scf_energy_list)), scf_energy_list))
        if len(res) < 2:
            return None
        return res

    def read_bandstructure(self):
        try:
            f = open(self.project_directory + self._working_dirctory + '/bands.out', 'r')
        except IOError:
            return None
        text = f.read().replace('-',' -')
        f.close()

        k_points = []
        energy_values = []
        found_line = False
        for line in text.split('\n'):
            line = line.strip()
            if line.strip().startswith('k ='):
                line_list = line.split()
                k_points.append([float(line_list[2]),float(line_list[3]),float(line_list[4])] )
                found_line = True
                continue
            if found_line and len(line)>0:
                e_split = line.split()
                e_numbers = [float(x) for x in e_split]
                energy_values.append(e_numbers)
                found_line = False

        n_bands = len(energy_values[0])
        n_k_points = len(k_points)
        bands = []
        k_array = np.zeros(n_k_points)
        for i in range(1,n_k_points):
            k_array[i] = np.linalg.norm(np.array(k_points[i])-np.array(k_points[i-1])) + k_array[i-1]

        for i in range(n_bands):
            band = np.zeros((n_k_points,2))
            band[:,0] = k_array

            e_band = [x[i] for x in energy_values]
            band[:,1] = e_band
            bands.append(band)



        return sst.BandStructure(bands,1,1)



    def read_gw_bandstructure(self, filename='BAND-QP.OUT'):
        raise NotImplementedError

    def read_phonon_bandstructure(self):
        raise NotImplementedError

    def read_optical_spectrum(self):
        raise NotImplementedError

    def read_ks_state(self):
        raise NotImplementedError

    def calculate_ks_density(self, crystal_structure, bs_point, grid='40 40 40'):
        raise NotImplementedError

    def kill_engine(self):
        try:
            self.engine_process.kill()
        except Exception as e:
            print(e)

    def is_engine_running(self, tasks=None):
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
        return True

    def _make_input_file(self,filename='scf.in'):
        if not os.path.isdir(self.project_directory + self._working_dirctory):
            os.mkdir(self.project_directory + self._working_dirctory)
        f = open(self.project_directory+self._working_dirctory+'/'+filename,'w')
        return f

    def _add_scf_to_file(self,file,crystal_structure,calculation='scf',band_points=None):
        if calculation == 'bands' and band_points is None:
            raise Exception('If calculation is bands you need to supply band points')

        control_options = {x: self.scf_options[x] for x in ['prefix']}
        control_options['verbosity'] = 'high'
        control_options['pseudo_dir'] = self.project_directory + self.pseudo_directory
        control_options['outdir'] = self.project_directory + self._working_dirctory
        control_options['calculation'] = calculation
        self._write_block(file, '&control',control_options )
        system_options = {}
        system_options['ecutwfc'] = float(self.scf_options['ecutwfc'])
        system_options['ibrav'] = 0
        system_options['nat'] = crystal_structure.n_atoms
        system_options['nbnd'] = int(self.scf_options['nbnd'])
        n_typ = set(crystal_structure.atoms[:,3])
        system_options['ntyp'] = len(n_typ)
        self._write_block(file,'&system',system_options)
        self._write_block(file,'&electrons',{})
        file.write('ATOMIC_SPECIES\n')
        for specie in n_typ:
            file.write(p_table[specie] + " {0:1.5f}".format(atomic_mass[specie]) +' '+ p_table[specie]+'.pseudo\n')

        file.write('ATOMIC_POSITIONS crystal\n')
        crystal_structure.atoms = crystal_structure.atoms[np.argsort(crystal_structure.atoms[:,3]),:]
        for i in range(crystal_structure.n_atoms):
            atom = crystal_structure.atoms[i,:]
            coords = atom[:3]
            specie = p_table[atom[3]]
            file.write(specie+' {0:1.5f} {1:1.5f} {2:1.5f}\n'.format(*coords))

        file.write('CELL_PARAMETERS bohr\n')
        for i in range(3):
            unit_vec = crystal_structure.lattice_vectors[i,:]
            file.write('{0:1.5f} {1:1.5f} {2:1.5f}\n'.format(*unit_vec))
        if calculation == 'scf':
            file.write('K_POINTS (automatic) \n')
            file.write(self.scf_options['k points']+ ' '+self.scf_options['k point shift'])
        elif calculation == 'bands':
            file.write('K_POINTS {crystal_b} \n')
            file.write('  '+str(len(band_points))+'\n')
            for band_point,label in band_points:
                file.write(' {0:1.5f} {1:1.5f} {2:1.5f} '.format(*band_point)+self.scf_options['k points band']+' !'+label+'\n')


    def _start_engine(self,filename='scf.in'):
        os.chdir(self.project_directory + self._working_dirctory)
        if self.custom_command_active:
            command = ['bash', self.custom_command]
        else:
            command = self._engine_command


        outname = filename.split('.')[0] + '.out'
        final_command = [None]
        final_command[0] = command[0] + ' <'+filename+' >'+outname

        self.engine_process = subprocess.Popen(final_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=True)
        os.chdir(self.project_directory)

    def _is_engine_running_custom_command(self,tasks):
        raise NotImplementedError

    def _write_block(self,file,block_name,options):
        file.write(block_name+'\n')
        for key,value in options.items():
            if type(value) == int:
                file.write('   '+key+'='+str(value)+'\n')
            elif type(value) == float:
                file.write('   '+key+'='+'{0:1.5f}'.format(value)+'\n')
            elif type(value) == str or isinstance(value, string_types):
                file.write('   '+key + '=' +"'"+value+"'" + '\n')
            else:
                raise Exception('bad type for option')

        file.write('/\n')

    def _correct_types(self):
        for key,value in self.scf_options_non_string_type.items():
            self.scf_options[key] = value(self.scf_options[key])

if __name__ == '__main__':
    atoms = np.array([[0, 0, 0, 6], [0.25, 0.25, 0.25, 6]])
    unit_cell = 6.719 * np.array([[0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    crystal_structure = sst.CrystalStructure(unit_cell, atoms)

    handler = Handler()
    handler.project_directory = "/home/jannick/OpenDFT_projects/test_qe"
    handler.scf_options['ecutwfc'] = 20.0
    handler.start_ground_state(crystal_structure,band_structure_points=((np.array([0,0,0]),'gamma'),(np.array([0.5,0.5,0.5]),'W')))
    handler.read_bandstructure()