import numpy as np
import solid_state_tools as sst
import periodictable as pt
import subprocess
import os
import time
import re
import threading
from six import string_types
from shutil import copyfile

atomic_mass = pt.mass
p_table = {i: el.__repr__() for i, el in enumerate(pt.elements)}
p_table_rev = {el.__repr__(): i for i, el in enumerate(pt.elements)}
atomic_mass = {i: el.mass for i, el in enumerate(pt.elements)}

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
    def __init__(self):
        self.engine_name = 'abinit'
        self.default_extension = '.xml'
        self._engine_command = ["abinit"]
        self.working_dirctory = '/abinit_files/'
        self.pseudo_directory = '/pseudos/'
        self.engine_process = None
        self.info_file = 'input.log'
        self.info_text = """ABINIT is a package whose main program allows one to find the total energy, charge density and electronic structure of systems made of electrons and nuclei (molecules and periodic solids) within Density Functional Theory (DFT), using pseudopotentials (or PAW atomic data) and a planewave basis. 
        ABINIT also optimize the geometry according to the DFT forces and stresses, or perform molecular dynamics simulations using these forces, 
        or generate phonons, Born effective charges, and dielectric tensors, based on Density-Functional Perturbation Theory, and many more properties. 
        Excited states and spectra can be computed within the Many-Body Perturbation Theory (the GW approximation and the Bethe-Salpeter equation). 
        DFT+U and Dynamical mean-field theory are available for strongly correlated materials. In addition to the main ABINIT code, different utility programs are provided.
ABINIT keywords are : capabilities, reliability, portability, documentation, with a "free software license" (short presentation of the ABINIT project - pdf document)."""

        self._filenames_tasks = {}
        self._timestamp_tasks = {}

        self._engine_version = self._get_engine_version()

        self.supported_methods = sst.ComputationalMethods(['periodic', 'scf', 'bandstructure'])

        self.project_directory = None
        self.input_filename = 'scf.in'

        self.current_input_file = self.input_filename
        self.current_output_file = self.info_file

        self.custom_command = ''
        self.custom_command_active = False
        self.dft_installation_folder = self.find_engine_folder()
        self.scf_options = {'ecut':'30.0','nstep':'30','toldfe1':'0.000001','ngkpt1':'5 5 5','nband':'10','ndivk':'30'}

        self.scf_options_tooltip = {}

        self.general_options = {'title': 'title'}
        self.bs_options = {}
        self.relax_options = {}
        self.relax_options_tooltip = {}

        self.gw_options = {}
        self.gw_options_tooltip = {}
        self.phonons_options = {}
        self.phonons_options_tooltip = {}
        self.optical_spectrum_options = {}
        self.optical_spectrum_options_tooltip = {}

        self.relax_file_timestamp = None

    def find_engine_folder(self):
        p = subprocess.Popen([search_command, 'abinit'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=shell_bool)
        res, err = p.communicate()
        res = res.decode()
        return res.strip()

    def parse_input_file(self, filename):
        raise NotImplementedError()

    def start_ground_state(self, crystal_structure, band_structure_points=None):
        pseudos = self._copy_default_pseudos(crystal_structure)
        self._make_files_file(pseudos)
        file = self._make_input_file()
        self._add_scf_to_file(file, crystal_structure,band_points=band_structure_points)
        file.close()
        self._start_engine()


    def start_optical_spectrum(self, crystal_structure):
        raise NotImplementedError

    def start_gw(self, crystal_structure, band_structure_points=None):
        raise NotImplementedError

    def start_phonon(self, crystal_structure, band_structure_points):
        raise NotImplementedError

    def start_relax(self, crystal_structure):
        if crystal_structure.n_atoms // 2 >= self.scf_options['nbnd']:
            raise Exception('Too few bands')
        self._copy_default_pseudos(crystal_structure)
        file = self._make_input_file()
        self._add_scf_to_file(file, crystal_structure, calculation='relax')
        file.close()
        self._start_engine()

    def load_relax_structure(self):
        file = self.project_directory + self.working_dirctory + self.info_file
        if not os.path.isfile(file):
            return None
        if self.relax_file_timestamp is not None and os.path.getmtime(file) == self.relax_file_timestamp:
            return None
        self.relax_file_timestamp = os.path.getmtime(file)

        with open(file) as f:
            text = f.readlines()

        matched_lines = [i for i, line in enumerate(text) if 'atomic_positions' in line.lower()]
        if len(matched_lines) == 0:
            return None
        highest_index = max(matched_lines)
        species = []
        coords = []
        for line in text[highest_index + 1:]:
            try:
                res = line.split()
                species.append(p_table_rev[res[0]])
                coords.append(np.array([float(x) for x in res[1:]]))
            except Exception:
                if len(res) == 0:
                    continue
                else:
                    break

        atoms = np.zeros((len(species), 4))
        for i, atom in enumerate(zip(species, coords)):
            atoms[i, :3] = atom[1]
            atoms[i, 3] = atom[0]

        matched_line = [i for i, line in enumerate(text) if 'lattice parameter' in line.lower()]
        line = text[matched_line[0]]
        res = line.split('=')[1].replace('a.u.', '')
        a = float(res)

        matched_line = [i for i, line in enumerate(text) if 'a(1) =' in line.lower()][0]
        lattice_vectors = np.zeros((3, 3))

        for i in range(matched_line, matched_line + 3):
            line = text[i]
            res = line.split('=')[1]
            res = res.replace('(', '').replace(')', '').split()
            res = np.array([float(x) for x in res])
            lattice_vectors[i - matched_line, :] = res * a

        return sst.CrystalStructure(lattice_vectors, atoms)

    def read_scf_status(self):
        try:
            f = open(self.project_directory + self.working_dirctory + self.info_file, 'r')
        except IOError:
            return None
        info_text = f.read()
        f.close()
        iteration_number = []
        scf_energy_list = []
        matches = re.findall(r"ETOT[\s\t]*\d+[\s\t]*[-+]?\d*\.\d+", info_text)
        for match in matches:
            ms = match.split()
            scf_energy_list.append(float(ms[2]))
            iteration_number.append(int(ms[1]))
        res = np.array(zip(iteration_number, scf_energy_list))
        if len(res) < 2:
            return None
        return res

    def read_bandstructure(self, special_k_points=None):
        try:
            f = open(self.project_directory + self.working_dirctory + '/scf_xo_DS2_EIG', 'r')
        except IOError:
            return None
        text = f.read()
        f.close()

        k_points = []
        energy_values = []
        found_line = False
        special_k_point_initial = []
        e_numbers = []
        for line in text.split('\n'):
            line = line.strip()
            if line.strip().startswith('kpt#'):
                if len(e_numbers)>0:
                    energy_values.append(e_numbers)

                e_numbers = []
                line_list = line.split()
                read_k_point = [float(line_list[7]), float(line_list[8]), float(line_list[9])]
                k_points.append(read_k_point)

                if special_k_points is not None:
                    for k_point, label in special_k_points:
                        if np.linalg.norm(np.array(read_k_point) - k_point) < 0.001:
                            special_k_point_initial.append([len(k_points) - 1, label])
                            break
                found_line = True
                continue
            if len(line) > 0 and found_line:
                e_split = line.split()
                e_numbers.extend([float(x) for x in e_split])

        if len(e_numbers) > 0:
            energy_values.append(e_numbers)

        n_bands = len(energy_values[0])
        n_k_points = len(k_points)
        bands = []
        k_array = np.zeros(n_k_points)
        for i in range(1, n_k_points):
            k_array[i] = np.linalg.norm(np.array(k_points[i]) - np.array(k_points[i - 1])) + k_array[i - 1]

        for i in range(n_bands):
            band = np.zeros((n_k_points, 2))
            band[:, 0] = k_array

            e_band = [x[i] for x in energy_values]
            band[:, 1] = e_band
            bands.append(band)

        special_k_points_out = [[k_array[i], label] for i, label in special_k_point_initial]

        try:
            f = open(self.project_directory + self.working_dirctory + self.info_file, 'r')
            info_text = f.read()
            f.close()

            matches = re.findall(r"nelect[\s\t]*=[\s\t]*[-+]?\d*\.\d+", info_text)
            match = matches[0]

            n_electrons = int(float(match.split('=')[1]))

            valence_bands = [band for i,band in enumerate(bands) if i<n_electrons//2]
            cond_bands = [band for i, band in enumerate(bands) if i >= n_electrons // 2]

            evalence = max( [band[:,1].max() for band in valence_bands] )
            econd = min([band[:,1].min() for band in cond_bands])

            efermi = evalence + (econd-evalence)/2
            for band in bands:
                band[:, 1] = band[:, 1] - efermi

        except IOError:
            pass



        return sst.BandStructure(bands, special_k_points=special_k_points_out)

    def read_gw_bandstructure(self, filename='BAND-QP.OUT'):
        raise NotImplementedError

    def read_phonon_bandstructure(self):
        raise NotImplementedError

    def read_optical_spectrum(self):
        raise NotImplementedError

    def read_ks_state(self):
        raise NotImplementedError


    def calculate_ks_density(self, crystal_structure, bs_point):
        raise NotImplementedError

    def calculate_electron_density(self, crystal_structure):
        raise NotImplementedError

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

    def reset_to_defaults(self):
        default_handler = Handler()
        self.scf_options.update(default_handler.scf_options)
        self.gw_options.update(default_handler.gw_options)
        self.optical_spectrum_options.update(default_handler.optical_spectrum_options)
        self.relax_options.update(default_handler.relax_options)
        self.phonons_options.update(default_handler.phonons_options)

    def _make_input_file(self, filename='scf.in'):
        if not os.path.isdir(self.project_directory + self.working_dirctory):
            os.mkdir(self.project_directory + self.working_dirctory)
        f = open(self.project_directory + self.working_dirctory + '/' + filename, 'w')
        return f

    def _add_scf_to_file(self, file, crystal_structure, band_points=None):

        if band_points is not None:
            file.write('ndtset 2\n')
        else:
            file.write('ndtset 1\n')


        file.write('# Definition of the unit cell\n')
        file.write('acell 1.0 1.0 1.0\n')
        file.write('rprim {0:1.10f} {1:1.10f} {2:1.10f}   {3:1.10f} {4:1.10f} {5:1.10f}   {6:1.10f} {7:1.10f} {8:1.10f}\n\n'.format(*np.reshape(crystal_structure.lattice_vectors,9)))

        atom_species = set(crystal_structure.atoms[:,3])
        n_atom = crystal_structure.atoms.shape[0]

        file.write('# Definition of atoms\n')
        file.write('ntypat {0:d}\n'.format(len(atom_species)))
        file.write('znucl')
        for atom_specie in atom_species:
            file.write(' {0:d}'.format(int(atom_specie)))
        file.write('\n')

        file.write('natom {0:d}\n'.format(n_atom))
        file.write('typat ')
        for i in range(n_atom):
            file.write('{0:d} '.format(list(atom_species).index(crystal_structure.atoms[i,3])+1))
        file.write('\n')
        file.write('xred\n')
        for i in range(n_atom):
            file.write('{0:1.10f} {1:1.10f} {2:1.10f}\n'.format(*crystal_structure.atoms[i,:3]))

        file.write('\n# Definition of the k grid\n')
        file.write('kptopt1 1\n')
        file.write('nshiftk1 4\n')
        file.write("""shiftk1 0.5 0.5 0.5
       0.5 0.0 0.0
       0.0 0.5 0.0
       0.0 0.0 0.5\n""")
        file.write('\n# All other options\n')
        for key,value in self.scf_options.items():
            if key not in ['nband','ndivk']:
                file.write(key+' '+value+'\n')

        file.write('prtden1 1\n')

        if band_points is not None:
            file.write("""\n#Dataset 2 : the band structure
iscf2    -2
getden2  -1
""")
            file.write('kptopt2 -{0:d}\n'.format(len(band_points)-1))
            file.write('nband2 ' + self.scf_options['nband']+'\n')
            file.write('ndivk2 ' + (len(band_points)-1) *(self.scf_options['ndivk']+' '))
            file.write('\nkptbounds2 ')
            for band_point in band_points:
                file.write('    {0:1.10f} {1:1.10f} {2:1.10f}\n'.format(*band_point[0]))
            file.write("""tolwfr2  1.0d-12\nenunit2  1   """)

    def _start_engine(self, filename='input.files'):
        os.chdir(self.project_directory + self.working_dirctory)
        if self.custom_command_active:
            command = ['bash', self.custom_command]
        else:
            command = self._engine_command

        outname = filename.split('.')[0] + '.log'
        final_command = [None]
        final_command[0] = command[0] + ' <' + filename + ' >' + outname

        self.engine_process = subprocess.Popen("exec " + final_command[0], stdout=subprocess.PIPE,
                                               stderr=subprocess.PIPE, shell=True)
        os.chdir(self.project_directory)

    def _is_engine_running_custom_command(self, tasks):
        raise NotImplementedError

    def _get_engine_version(self):
        p = subprocess.Popen(['abinit','--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=shell_bool)
        res, err = p.communicate()
        res = res.decode()
        version_list = res.strip().split('.')
        version_int_list = [int(x) for x in version_list]
        return version_int_list


    def _copy_default_pseudos(self, crystal_structure):
        atoms = set(crystal_structure.atoms[:,3])
        atoms_names = [p_table[atom] for atom in atoms]
        installation_folder = os.path.dirname(__file__)

        if not os.path.isdir(self.project_directory+self.pseudo_directory):
            os.mkdir(self.project_directory+self.pseudo_directory)

        pseudo_files = []
        for atom in atoms_names:
            file = atom+'.psp8'
            pseudo_files.append(file)
            filepath = self.project_directory+self.pseudo_directory+file
            if not os.path.isfile(filepath):
                copyfile(installation_folder+'/data/pseudos/abinit/'+file,filepath)

            if self._engine_version[0]<7 or (self._engine_version[0]==7 and self._engine_version[1]<10):
                with open(filepath, 'r') as f:
                    filedata = f.readlines()
                filedata[5] = filedata[5].lstrip()
                filedata[5] = '0' + filedata[5][1:]

                with open(filepath, 'w') as f:
                    f.writelines(filedata)

        return pseudo_files


    def _make_files_file(self,pseudos):
        if not os.path.isdir(self.project_directory + self.working_dirctory):
            os.mkdir(self.project_directory + self.working_dirctory)
        with open(self.project_directory+self.working_dirctory+'/input.files','w') as f:
            f.write('scf.in\n')
            f.write('scf.out\n')
            f.write('scf_xi\n')
            f.write('scf_xo \n')
            f.write('scf_x\n')
            for pseudo in pseudos:
                f.write('../pseudos/'+pseudo+'\n')


if __name__ == '__main__':
    atoms = np.array([[0, 0, 0, 14], [0.25, 0.25, 0.25, 14]])
    unit_cell = np.array([[3.3, 0.0, 3.3], [3.3, 3.3, 0.0], [0, 3.3, 3.3]])
    crystal_structure = sst.CrystalStructure(unit_cell, atoms,scale=6.6)

    handler = Handler()
    handler.project_directory = "/home/jannick/OpenDFT_projects/test_abinit"
    # handler.scf_options['ecutwfc'] = 20.0
    band_structure_points = ((np.array([0, 0, 0]), 'gamma'), (np.array([0.5, 0.5, 0.5]), 'W'), (np.array([0.0, 0.0, 0.5]), 'Z'), (np.array([0.5, 0.0, 0.0]), 'X'), (np.array([0.0, 0.5, 0.0]), 'Y'))

    handler._copy_default_pseudos(crystal_structure)

    # handler.start_ground_state(crystal_structure, band_structure_points=band_structure_points)
    # while handler.is_engine_running():
    #     time.sleep(0.3)
    # res = handler.read_scf_status()

    band_structure = handler.read_bandstructure(special_k_points=band_structure_points)
    version = handler._get_engine_version()