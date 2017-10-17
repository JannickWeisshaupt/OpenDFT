import numpy as np
import solid_state_tools as sst
import periodictable as pt
import subprocess
import os
import re

atomic_mass = pt.mass
p_table = {i: el.__repr__() for i, el in enumerate(pt.elements)}
p_table_rev = {el.__repr__(): i for i, el in enumerate(pt.elements)}
atomic_mass = {i: el.mass for i, el in enumerate(pt.elements)}

hartree = 27.211
bohr = 0.52917721

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
        self.engine_name = 'nwchem'
        self.default_extension = ''
        self._engine_command = ["nwchem"]
        self.working_dirctory = '/nwchem_files/'
        self.pseudo_directory = None
        self.engine_process = None
        self.info_file = 'scf.out'
        self.info_text = """NWChem aims to provide its users with computational chemistry tools that are scalable both in their ability to treat large scientific computational chemistry problems efficiently, and in their use of available parallel computing resources from high-performance parallel supercomputers to conventional workstation clusters.

NWChem software can handle:

 - Biomolecules, nanostructures, and solid-state
 - From quantum to classical, and all combinations
 - Ground and excited-states
 - Gaussian basis functions or plane-waves
 - Scaling from one to thousands of processors
 - Properties and relativistic effects 

NWChem is actively developed by a consortium of developers and maintained by the EMSL located at the Pacific Northwest National Laboratory (PNNL) in Washington State. 
Researchers interested in contributing to NWChem should review the Developers page. The code is distributed as open-source under the terms of the Educational Community License version 2.0 (ECL 2.0).

The NWChem development strategy is focused on providing new and essential scientific capabilities to its users in the areas of kinetics and dynamics of chemical transformations, chemistry at interfaces and in the condensed phase, and enabling innovative and integrated research at EMSL. 
At the same time continued development is needed to enable NWChem to effectively utilize architectures of tens of petaflops and beyond. """

        self._filenames_tasks = {}
        self._timestamp_tasks = {}

        self.supported_methods = sst.ComputationalMethods(
            ['non-periodic', 'scf','relax'])

        self.project_directory = None
        self.input_filename = 'scf.in'
        self.custom_command = ''
        self.custom_command_active = False
        self.dft_installation_folder = self.find_engine_folder()
        self.scf_options = {'basis': 'cc-pvdz','method':'scf','spin state':'singlet','spin restriction':'RHF'}
        self.scf_options_tooltip = {}

        self.general_options = {'title': 'title'}
        self.bs_options = {}
        self.relax_options = {}
        self.gw_options = {}
        self.gw_options_tooltip = {}
        self.phonons_options = {}
        self.phonons_options_tooltip = {}
        self.optical_spectrum_options = {}

        self.relax_file_timestamp = None

    def find_engine_folder(self):
        p = subprocess.Popen([search_command, 'nwchem'], stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=shell_bool)
        res, err = p.communicate()
        res = res.decode()
        res = res.split('bin')[0]
        return res.strip()

    def parse_input_file(self, filename):
        raise NotImplementedError()

    def start_ground_state(self, crystal_structure, band_structure_points=None):
        file = self._make_input_file()
        self._add_scf_to_file(file,crystal_structure)
        file.close()
        self._start_engine()

        # if band_structure_points is not None:
        #     def run_bs():
        #         while self.is_engine_running():
        #             time.sleep(0.001)
        #         file = self._make_input_file(filename='bands.in')
        #         self._add_scf_to_file(file,crystal_structure,calculation='bands',band_points=band_structure_points)
        #         file.close()
        #         self._start_engine(filename='bands.in')
        #
        #
        #     t = threading.Thread(target=run_bs)
        #     t.start()

    def start_optical_spectrum(self, crystal_structure):
        raise NotImplementedError

    def start_gw(self, crystal_structure, band_structure_points=None):
        raise NotImplementedError

    def start_phonon(self, crystal_structure, band_structure_points):
        raise NotImplementedError

    def start_relax(self, crystal_structure):
        file = self._make_input_file()
        self._add_scf_to_file(file,crystal_structure,calculation='optimize')
        file.close()
        self._start_engine()

    def load_relax_structure(self):
        file = self.project_directory+self.working_dirctory+self.info_file
        if not os.path.isfile(file):
            return None
        if self.relax_file_timestamp is not None and os.path.getmtime(file) == self.relax_file_timestamp:
            return None
        self.relax_file_timestamp = os.path.getmtime(file)

        with open(file) as f:
            text = f.readlines()

        matched_lines = [i for i,line in enumerate(text) if 'output coordinates in a.u.' in line.lower()]
        if len(matched_lines) == 0:
            return None
        highest_index = max(matched_lines)
        species = []
        coords = []
        for line in text[highest_index+4:]:
            try:
                res = line.split()
                species.append(p_table_rev[res[1].title()])
                coords.append(np.array([float(x) for x in res[3:]]))
            except Exception:
                if len(res)==0:
                    continue
                else:
                    break

        atoms = np.zeros((len(species),4))
        for i,atom in enumerate(zip(species,coords)):
            atoms[i,:3] = atom[1]
            atoms[i,3] =atom[0]

        return sst.MolecularStructure(atoms)

    def read_scf_status(self):
        try:
            f = open(self.project_directory + self.working_dirctory + self.info_file, 'r')
        except IOError:
            return None
        info_text = f.read()
        f.close()

        scf_energy_list = []

        hits = [i for i,l in enumerate(info_text.split('\n')) if l.strip().startswith('iter')]
        if len(hits)==0:
            return None
        else:
            hit = hits[0]
        for line in info_text.split('\n')[hit+2:]:
            sline = line.split()
            if len(sline)==0:
                break
            scf_energy_list.append(float(sline[1]))

        matches = re.findall(r"Total SCF energy[\s\t]*=[\s\t]*[-+]?\d*\.\d+", info_text)
        for match in matches:
            ms = match.split('=')
            scf_energy_list.append(float(ms[1]))

        res = np.array(zip(range(1,len(scf_energy_list)+1), scf_energy_list))

        if len(res) < 2:
            return None
        return res

    def read_bandstructure(self,special_k_points=None):
        raise NotImplementedError

    def read_energy_diagram(self):
        try:
            f = open(self.project_directory + self.working_dirctory + self.info_file, 'r')
        except IOError:
            return None
        text = f.readlines()
        f.close()

        # hits = [i for i,l in enumerate(text) if 'final eigenvalues' in l.lower()]
        # if len(hits) == 0:
        #     return None
        # hit = hits[-1]
        #
        # energies = []
        # started_bool = False
        # for line in text[hit+1:]:
        #     if '-----' in line:
        #         continue
        #     try:
        #         sline = line.split()
        #         if len(sline)==0:
        #             continue
        #         energies.append(float(sline[1])*hartree)
        #         started_bool = True
        #
        #     except Exception:
        #         if started_bool:
        #             break
        #
        # labels = [str(i+1) for i,energie in enumerate(energies)]
        # hits = [i for i,l in enumerate(text) if 'final eigenvalues' in l.lower()]
        # if len(hits) == 0:
        #     return None
        # hit = hits[-1]
        #
        # energies = []
        # started_bool = False
        # for line in text[hit+1:]:
        #     if '-----' in line:
        #         continue
        #     try:
        #         sline = line.split()
        #         if len(sline)==0:
        #             continue
        #         energies.append(float(sline[1])*hartree)
        #         started_bool = True
        #
        #     except Exception:
        #         if started_bool:
        #             break
        #
        # labels = [str(i+1) for i,energie in enumerate(energies)]

        hits = [i for i,l in enumerate(text) if 'final molecular orbital analysis' in l.lower()]
        if len(hits) == 0:
            return None
        hit = hits[-1]

        energies = []
        occupations = []
        for line in text[hit+1:]:
            try:
                if not line.strip().startswith('Vector'):
                    continue
                sline = line.split('E=')
                if len(sline)==0:
                    continue
                energies.append(float(sline[1].strip().replace('D','e'))*hartree)
                sline2 = line.split('Occ=')[1]
                sline3 = sline2.split('E=')[0]
                sline3 = sline3.replace('D','e')
                occupations.append(float(sline3))
            except Exception:
                pass


        labels = ['' for i,energie in enumerate(energies)]
        return sst.EnergyDiagram(energies,labels,occupations=occupations)

    def read_gw_bandstructure(self, filename='BAND-QP.OUT'):
        raise NotImplementedError

    def read_phonon_bandstructure(self):
        raise NotImplementedError

    def read_optical_spectrum(self):
        raise NotImplementedError

    def read_ks_state(self):

        with open(self.project_directory+self.working_dirctory+ '/chargedensity.cube') as f:
            text = f.readlines()
        origin = np.array(text[2].split()[1:],dtype=np.float)
        lattice_vecs = np.zeros((3,3))
        for i in range(3):
            line = text[3+i]
            coords = np.array(line.split()[1:],dtype=np.float)
            steps = int(line.split()[0])
            lattice_vecs[i,:3] = coords*(steps-1)

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

        if data.min() == data.max():
            return None

        data = data / data.max()
        return sst.MolecularDensity(data,lattice_vecs,origin)

    def calculate_ks_density(self, crystal_structure, bs_point):
        file = self._make_input_file()
        file.write('title '+'"'+self.general_options['title']+'"\n')
        self._add_geometry(file,crystal_structure,auto=False)
        self._add_basis(file,crystal_structure)
        # self._add_scf_field_to_file(file)
        self._add_dplot_to_file(file,crystal_structure,orbital=bs_point[1])
        file.write('task dplot\n')
        file.close()
        self._start_engine()

    def calculate_electron_density(self,crystal_structure):
        file = self._make_input_file()
        file.write('title '+'"'+self.general_options['title']+'"\n')
        self._add_geometry(file,crystal_structure,auto=False)
        self._add_basis(file,crystal_structure)
        # self._add_scf_field_to_file(file)
        self._add_dplot_to_file(file,crystal_structure)
        file.write('task dplot\n')
        file.close()
        self._start_engine()

    def kill_engine(self):
        try:
            self.engine_process.kill()
            # os.killpg(os.getpgid(self.engine_process.pid), signal.SIGTERM)
        except Exception as e:
            print(e)

    def will_scf_run(self):
        return True

    def reset_to_defaults(self):
        default_handler = Handler()
        self.scf_options.update(default_handler.scf_options)
        self.gw_options.update(default_handler.gw_options)
        self.optical_spectrum_options.update(default_handler.optical_spectrum_options)
        self.relax_options.update(default_handler.relax_options)
        self.phonons_options.update(default_handler.phonons_options)

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

    def _add_scf_to_file(self, file, crystal_structure, calculation='scf', band_points=None):
        file.write('title '+'"'+self.general_options['title']+'"\n')
        self._add_geometry(file,crystal_structure)
        self._add_basis(file,crystal_structure)
        self._add_scf_field_to_file(file)
        file.write('task ' + self.scf_options['method'])
        if calculation == 'optimize':
            file.write(' '+calculation)
        file.write('\n')

    def _make_input_file(self, filename='scf.in'):
        if not os.path.isdir(self.project_directory + self.working_dirctory):
            os.mkdir(self.project_directory + self.working_dirctory)
        f = open(self.project_directory + self.working_dirctory + '/' + filename, 'w')
        return f

    def _start_engine(self, filename='scf.in'):
        os.chdir(self.project_directory + self.working_dirctory)
        if self.custom_command_active:
            command = ['bash', self.custom_command]
        else:
            command = self._engine_command

        outname = filename.split('.')[0] + '.out'
        final_command = [None]
        final_command[0] = command[0] +' '+ filename + ' >' + outname

        self.engine_process = subprocess.Popen("exec " + final_command[0], stdout=subprocess.PIPE,
                                               stderr=subprocess.PIPE, shell=True)
        os.chdir(self.project_directory)

    def _add_geometry(self,file,crystal_structure,auto=False):
        if not auto:
            auto_string = 'noautoz nocenter noautosym'
        else:
            auto_string = ''
        file.write('geometry units atomic '+auto_string+'\n')

        for atom in crystal_structure.atoms:
            file.write('    '+p_table[atom[3]].lower())
            file.write(" {0:1.9f} {1:1.9f} {2:1.9f}\n".format(*atom[:3]))
        file.write('end\n')

    def _add_basis(self,file,crystal_structure):
        file.write('basis\n')
        file.write('    *'+ ' library ' + self.scf_options['basis']+'\n')
        # for atom in crystal_structure.atoms:
        #     file.write('    '+p_table[atom[3]].lower() + ' library ' + self.scf_options['basis']+'\n')
        file.write('end\n')

    def _add_scf_field_to_file(self,file,input=False):
        file.write('scf\n')
        if input:
            file.write('    vectors input scf.movecs\n')
        file.write('   '+self.scf_options['spin state']+'\n')
        file.write('   '+self.scf_options['spin restriction']+'\n')

        file.write('end\n')

    def _add_dplot_to_file(self,file,crystal_structure,limits=None,orbital=None):
        if limits is None:
            edge = 0.5
            limits = np.zeros((3,3))
            coords = crystal_structure.calc_absolute_coordinates()[:,:3]
            for i in range(3):
                coord = coords[:,i]
                tr_range = coord.max()-coord.min()
                if tr_range < 0.5:
                    tr_range = 5
                limits[i,:] = [coord.min()-tr_range*edge,coord.max()+tr_range*edge,100]





        file.write('dplot\n')
        file.write('    TITLE '+self.general_options['title']+'\n')
        file.write('    vectors scf.movecs\n')
        file.write('    LimitXYZ\n')
        for limit in limits:
            limit[:2] = limit[:2]*bohr
            file.write('    {0:1.9f} {1:1.9f} {2:1.0f}\n'.format(*limit))


        file.write("""    spin total
    gaussian
    output chargedensity.cube\n""")
        if orbital is not None:
            file.write('    orbitals density; 1; '+str(orbital)+'\n')

        file.write('end\n')

if __name__ == '__main__':
    atoms = np.array([[0, 0, -1.2, 7], [0, 0, 1.2, 7]])
    crystal_structure = sst.MolecularStructure(atoms)

    handler = Handler()
    handler.project_directory = "/home/jannick/OpenDFT_projects/nwchem_test"
    # handler.start_ground_state(crystal_structure)
    # handler.read_scf_status()
    energies = handler.read_energy_diagram()

    # handler.calculate_electron_density(crystal_structure)
    # while handler.is_engine_running():
    #     time.sleep(0.1)
    # dens = handler.read_ks_state()

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
    # handler.calculate_ks_density(None,[1,1])
    # while handler.is_engine_running():
    #     time.sleep(0.1)
    # ks = handler.read_ks_state()
    # out = handler.load_relax_structure()