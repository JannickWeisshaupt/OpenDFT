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
from collections import OrderedDict
import pandas as pd
import signal

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
        self.engine_name = 'quantum espresso'
        self.default_extension = '.xml'
        self._engine_command = ["pw.x"]
        self.working_dirctory = '/quantum_espresso_files/'
        self.pseudo_directory = '/pseudos/'
        self.engine_process = None
        self.info_file = 'scf.out'
        self._filenames_tasks = {}

        self._timestamp_tasks = {}

        self.project_directory = None
        self._input_filename = 'input.xml'
        self.custom_command = ''
        self.custom_command_active = False
        self.dft_installation_folder = self.find_engine_folder()
        self.scf_options = {'ecutwfc':'30.0','ecutrho':'300.0','input_dft':'PBE','k points':'6 6 6','k point shift':'1 1 1','k points band':'30','nbnd':'10',
                            'diagonalization':'david','conv_thr':'1e-8','mixing_mode':'plain','mixing_beta':'0.7'}


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
        p = subprocess.Popen(['which', 'pw.x'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        res, err = p.communicate()
        res = res.decode()
        res = res.split('bin')[0]
        return res.strip()

    def parse_input_file(self, filename):
        raise NotImplementedError()

    def start_ground_state(self, crystal_structure, band_structure_points=None):
        # self._correct_types()
        if crystal_structure.n_atoms//2 >= self.scf_options['nbnd']:
            raise Exception('Too few bands')

        file = self._make_input_file()
        self._add_scf_to_file(file,crystal_structure)
        file.close()
        self._start_engine()

        if band_structure_points is not None:
            def run_bs():
                while self.is_engine_running():
                    time.sleep(0.001)
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
            f = open(self.project_directory + self.working_dirctory + self.info_file, 'r')
        except IOError:
            return None
        info_text = f.read()
        f.close()


        scf_energy_list = []
        matches = re.findall(r"total energy[\s\t]*=[\s\t]*[-+]?\d*\.\d+", info_text)
        for match in matches:
            ms = match.split('=')
            scf_energy_list.append(float(ms[1]))

        res = np.array(zip(range(1,len(scf_energy_list)+1), scf_energy_list))
        if len(res) < 2:
            return None
        return res

    def read_bandstructure(self,special_k_points=None):
        try:
            f = open(self.project_directory + self.working_dirctory + '/bands.out', 'r')
        except IOError:
            return None
        text = f.read().replace('-',' -')
        f.close()

        k_points = []
        energy_values = []
        found_line = False
        special_k_point_initial = []
        for line in text.split('\n'):
            line = line.strip()
            if line.strip().startswith('k ='):
                line_list = line.split()
                read_k_point = [float(line_list[2]),float(line_list[3]),float(line_list[4])]
                k_points.append(read_k_point)

                if special_k_points is not None:
                    for k_point,label in special_k_points:
                        if np.linalg.norm( np.array(read_k_point) - k_point)<0.001:
                            special_k_point_initial.append([len(k_points)-1,label])
                            break

                found_line = True
                e_numbers = []
                continue
            if found_line and len(line)>0:
                e_split = line.split()
                e_numbers.extend([float(x) for x in e_split])
            elif found_line and len(line) == 0 and len(e_numbers)>0:
                energy_values.append(e_numbers)
                found_line = False

        matches = re.findall('number of electrons[\s\t]*=[\s\t]*[-+]?\d*\.\d+',text)
        n_electrons = int(float(matches[0].split('=')[1]))

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

        special_k_points_out = [[k_array[i],label] for i,label in special_k_point_initial]

        try:
            valence_bands = [band for i,band in enumerate(bands) if i<n_electrons//2]
            cond_bands = [band for i, band in enumerate(bands) if i >= n_electrons // 2]

            evalence = max( [band[:,1].max() for band in valence_bands] )
            econd = min([band[:,1].min() for band in cond_bands])

            efermi = evalence + (econd-evalence)/2
            for band in bands:
                band[:, 1] = band[:, 1] - efermi

        except Exception:
            pass

        return sst.BandStructure(bands,special_k_points=special_k_points_out)

    def read_gw_bandstructure(self, filename='BAND-QP.OUT'):
        raise NotImplementedError

    def read_phonon_bandstructure(self):
        raise NotImplementedError

    def read_optical_spectrum(self):
        raise NotImplementedError

    def read_ks_state(self):

        with open(self.project_directory+self.working_dirctory+ '/rho.dat') as f:
            text = f.readlines()
        total_res = []
        for line in text[2:]:
            res = line.split()
            if not '.' in res[0]:
                continue
            numbers = [float(x) for x in res]
            total_res.extend(numbers)

        n_grid = [int(text[i].split()[0]) for i in range(3,6)]

        l_data = np.array(total_res)
        data = l_data.reshape(n_grid,order='C')

        data = data / data.max()
        return sst.KohnShamDensity(data)

    def calculate_ks_density(self, crystal_structure, bs_point):
        f = self._make_input_file(filename='pp.in')
        inputpp_dic = {'prefix':self.general_options['title'],'outdir':self.project_directory + self.working_dirctory, 'plot_num':7, 'filplot': 'e_density', 'kpoint(1)':bs_point[0], 'kband(1)':bs_point[1]}
        self._write_block(f,'&inputpp',inputpp_dic)
        plot_dic = {'iflag':3,'output_format':6,'fileout':'rho.dat'}
        self._write_block(f,'&plot',plot_dic)
        f.close()
        self._start_pp_process()

    def calculate_electron_density(self,crystal_structure):
        f = self._make_input_file(filename='pp.in')
        inputpp_dic = {'prefix':self.general_options['title'],'outdir':self.project_directory + self.working_dirctory, 'plot_num':0, 'filplot': 'e_density'}
        self._write_block(f,'&inputpp',inputpp_dic)
        plot_dic = {'iflag':3,'output_format':6,'fileout':'rho.dat'}
        self._write_block(f,'&plot',plot_dic)
        f.close()
        self._start_pp_process()

    def _start_pp_process(self):
        command = 'exec pp.x<pp.in'
        os.chdir(self.project_directory + self.working_dirctory)
        self.engine_process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                               shell=True,preexec_fn=os.setpgrp)
        os.chdir(self.project_directory)

    def kill_engine(self):
        try:
            self.engine_process.kill()
            # os.killpg(os.getpgid(self.engine_process.pid), signal.SIGTERM)
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
        if not os.path.isdir(self.project_directory + self.working_dirctory):
            os.mkdir(self.project_directory + self.working_dirctory)
        f = open(self.project_directory + self.working_dirctory + '/' + filename, 'w')
        return f

    def _add_scf_to_file(self,file,crystal_structure,calculation='scf',band_points=None):
        if calculation == 'bands' and band_points is None:
            raise Exception('If calculation is bands you need to supply band points')

        control_options ={'prefix':self.general_options['title']}
        control_options['verbosity'] = 'high'
        control_options['pseudo_dir'] = self.project_directory + self.pseudo_directory
        control_options['outdir'] = self.project_directory + self.working_dirctory
        control_options['calculation'] = calculation
        self._write_block(file, '&control',control_options )
        system_options = {}
        system_options['ecutwfc'] = float(self.scf_options['ecutwfc'])
        system_options['ecutrho'] = float(self.scf_options['ecutrho'])
        system_options['ibrav'] = 0
        system_options['nat'] = crystal_structure.n_atoms
        system_options['nbnd'] = int(self.scf_options['nbnd'])

        electron_options = {'diagonalization':self.scf_options['diagonalization'],'conv_thr':float(self.scf_options['conv_thr']),'mixing_beta':float(self.scf_options['mixing_beta']),
                            'mixing_mode':self.scf_options['mixing_mode']}

        n_typ = set(crystal_structure.atoms[:,3])
        system_options['ntyp'] = len(n_typ)
        self._write_block(file,'&system',system_options)
        self._write_block(file,'&electrons',electron_options)
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
        os.chdir(self.project_directory + self.working_dirctory)
        if self.custom_command_active:
            command = ['bash', self.custom_command]
        else:
            command = self._engine_command


        outname = filename.split('.')[0] + '.out'
        final_command = [None]
        final_command[0] = command[0] + ' <'+filename+' >'+outname

        self.engine_process = subprocess.Popen("exec "+final_command[0], stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=True)
        os.chdir(self.project_directory)

    def _is_engine_running_custom_command(self,tasks):
        raise NotImplementedError

    def _write_block(self,file,block_name,options):
        file.write(block_name+'\n')
        for key,value in options.items():
            if type(value) == int:
                file.write('   '+key+'='+str(value)+'\n')
            elif type(value) == float:
                file.write('   '+key+'='+'{0:1.16f}'.format(value)+'\n')
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
    # handler.scf_options['ecutwfc'] = 20.0
    # band_structure_points = ((np.array([0,0,0]),'gamma'),(np.array([0.5,0.5,0.5]),'W'))
    # handler.start_ground_state(crystal_structure,band_structure_points=band_structure_points)
    #
    # coords = [x[0] for x in band_structure_points]
    # labels = [x[1] for x in band_structure_points]
    # new_coords = crystal_structure.convert_to_tpiba(coords)
    # band_structure_points_conv = zip(new_coords, labels)
    #
    # band_structure = handler.read_bandstructure(special_k_points=band_structure_points_conv)
    handler.calculate_ks_density(None,[1,1])
    while handler.is_engine_running():
        time.sleep(0.1)
    ks = handler.read_ks_state()