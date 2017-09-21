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
        self._engine_command = ["pw.x < scf.in >scf.out"]
        self._working_dirctory = '/quantum_espresso_files/'
        self._pseudo_dirctory = '/pseudos/'
        self.engine_process = None
        self._info_file = 'INFO.OUT'
        self._filenames_tasks = {}

        self._timestamp_tasks = {}

        self.project_directory = None
        self._input_filename = 'input.xml'
        self.custom_command = ''
        self.custom_command_active = False
        self.exciting_folder = self.find_engine_folder()
        self.scf_options = {'prefix':'title','ecutwfc':'12.0'}
        self.scf_options_non_string_type  = {'ecutwfc':float}


    def find_engine_folder(self):
        p = subprocess.Popen(['which', 'pw.x'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        res, err = p.communicate()
        res = res.decode()
        res = res.split('bin')[0]
        return res.strip()

    def parse_input_file(self, filename):
        raise NotImplementedError

    def start_ground_state(self, crystal_structure, band_structure_points=None):
        self._correct_types()
        file = self._make_input_file()
        self._add_scf_to_file(file,crystal_structure)
        file.close()
        self._start_engine()

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
        raise NotImplementedError

    def read_bandstructure(self):
        raise NotImplementedError

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

    def _make_input_file(self):
        if not os.path.isdir(self.project_directory + self._working_dirctory):
            os.mkdir(self.project_directory + self._working_dirctory)
        f = open(self.project_directory+self._working_dirctory+'/scf.in','w')
        return f

    def _add_scf_to_file(self,file,crystal_structure):
        control_options = {x: self.scf_options[x] for x in ['prefix']}
        control_options['pseudo_dir'] = self.project_directory + self._pseudo_dirctory
        control_options['outdir'] = self.project_directory + self._working_dirctory
        self._write_block(file, '&control',control_options )
        system_options = {x: self.scf_options[x] for x in ['ecutwfc']}
        system_options['ibrav'] = 0
        system_options['nat'] = crystal_structure.n_atoms
        n_typ = set(crystal_structure.atoms[:,3])
        system_options['ntyp'] = len(n_typ)
        self._write_block(file,'&system',system_options)
        self._write_block(file,'&electrons',{})
        file.write('ATOMIC_SPECIES\n')
        for specie in n_typ:
            file.write(p_table[specie] + " {0:1.5f}".format(atomic_mass[specie]) +' '+ p_table[specie]+'.pseudo\n')

        file.write('ATOMIC_POSITIONS\n')
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


    def _start_engine(self):
        os.chdir(self.project_directory + self._working_dirctory)
        if self.custom_command_active:
            command = ['bash', self.custom_command]
            # if tasks is not None:
            #     filenames = [self.filenames_tasks[task] for task in tasks if task != 'scf']
            #     for filename in filenames:
            #         os.remove(self.project_directory+self.working_dirctory+filename)
        else:
            command = self._engine_command

        self.engine_process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=True)
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
            elif type(value) == str:
                file.write('   '+key + '=' +"'"+value+"'" + '\n')
            else:
                raise Exception('bad type for option')

        file.write('/\n')

    def _correct_types(self):
        for key,value in self.scf_options_non_string_type:
            self.scf_options[key] = value(self.scf_options[key])

if __name__ == '__main__':
    atoms = np.array([[0, 0, 0, 6], [0.25, 0.25, 0.25, 6]])
    unit_cell = 6.719 * np.array([[0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    crystal_structure = sst.CrystalStructure(unit_cell, atoms)

    handler = Handler()
    handler.project_directory = "/home/jannick/OpenDFT_projects/test_qe"
    handler.start_ground_state(crystal_structure)