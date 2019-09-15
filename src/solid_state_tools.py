from __future__ import division
import numpy as np
import re
import periodictable as pt
from bisect import bisect
import time
from scipy.spatial import ConvexHull, Voronoi
import os
import copy

try:
    from little_helpers import find_data_file
except ImportError:
    from .little_helpers import find_data_file

from collections import OrderedDict
import itertools

from pymatgen.symmetry.bandstructure import HighSymmKpath,SpacegroupAnalyzer
from pymatgen.symmetry.analyzer import PointGroupAnalyzer
import pymatgen as mg
from pymatgen.ext.matproj import MPRester

API_KEY = "4m1ieH0XhKw2AzSL"

bohr = 0.52917721067
unit_mass = 1.660539040e-27

cov_radii = np.loadtxt(find_data_file('/data/cov_radii.dat')) / bohr

p_table = {i: el.__repr__() for i, el in enumerate(pt.elements)}
p_table_rev = {el.__repr__(): i for i, el in enumerate(pt.elements)}
elemental_masses = [pt.mass.mass(el) for el in pt.elements]

def remove_duplicates(data, treshold=0.01):
    if len(data) == 0:
        return np.array([])
    new_data = []
    n, m = data.shape
    for iii in range(n):

        tr_data = np.linalg.norm(data[iii + 1:, :] - data[iii, :], axis=1)
        if not any(tr_data < treshold):
            new_data.append(data[iii, :])
    return np.array(new_data)


class MolecularStructure(object):
    """Represents a molecular structure
    :param  atoms -- a (Nx4) numpy array with [x,y,z,type] as rows. Type is an integer with the atomic number (H=1,He=2 etc.).


    """
    def __init__(self, atoms, scale=1.0):
        self.atoms = np.array(atoms, dtype=np.float)  # np array with [x,y,z,type] type is number in periodic system
        self.atoms[:, :3] = self.atoms[:, :3] * scale
        self.n_atoms = atoms.shape[0]
        self.scale = scale  # This is just bonus info. Do not use this here. Only for editing

    def calc_absolute_coordinates(self, repeat=[1, 1, 1],edges=False):
        """Returns the cartesion coordinates"""
        return self.atoms

    def find_bonds(self, abs_coords):
        """Returns a list of bonds between atom i and j, i.e. ((i,j),...).
The bonds are calculated according to covalent radii from the literature.

        :return bonds:  a tuple of of two element tuples which represent the connections or bonds


"""
        n_atoms = abs_coords.shape[0]
        abs_coords_pure = abs_coords[:, :3]

        dist_mat = np.zeros((n_atoms, n_atoms))
        bonds = []
        for j1 in range(n_atoms):
            for j2 in range(j1 + 1, n_atoms):

                dist = np.linalg.norm(abs_coords_pure[j1, :] - abs_coords_pure[j2, :])
                if dist < (cov_radii[int(abs_coords[j1, 3])] + cov_radii[int(abs_coords[j2, 3])]) * 1.3:
                    dist_mat[j1, j2] = dist
                    bonds.append((j1, j2))
        return tuple(bonds)

    def symmetry_information(self):
        """
        Returns symmetry information of the molecule as for example its point group

        :return info: a dictionary with symmetry informations
        """
        mol = mg.Molecule(self.atoms[:, 3], self.atoms[:, :3])
        analyzer = PointGroupAnalyzer(mol)
        info = {'point group':analyzer.get_pointgroup()}
        return info


class CrystalStructure(object):
    """This class represents a crystal structure

    :param lattice_vectors: a 3x3 numpy array with the lattice vectors as rows
    :param atoms: a Nx4 numpy array the atomic positions in [:,:3] as either relative (crystal) or cartesion coordinates
    :param relative_coords:
    :param scale: A scale of the whole unit cell. This parameter is only used to nicely format the lattice vectors. The lattice_vectors argument must still contain the correct unit vectors.

    """
    def __init__(self, lattice_vectors, atoms, relative_coords=True, scale=1.0):
        self._lattice_vectors = np.array(lattice_vectors, dtype=np.float)  # tuple of np.arrays

        self.calculate_inv_lattice()

        self.atoms = np.array(atoms, dtype=np.float)  # np array with [x,y,z,type] type is number in periodic system
        self.n_atoms = atoms.shape[0]
        self.scale = scale  # This is just bonus info. Do not use this here. Only for editing

        if relative_coords:
            self.atoms[:, :3] = np.mod(self.atoms[:, :3], 1)
        else:
            inv_lattice_vecs = np.linalg.inv(self.lattice_vectors.T)
            for i in range(self.n_atoms):
                self.atoms[i, :3] = np.mod(np.dot(inv_lattice_vecs, self.atoms[i, :3].T),1)

    @property
    def lattice_vectors(self):
        return self._lattice_vectors

    @lattice_vectors.setter
    def lattice_vectors(self, value):
        self._lattice_vectors = value
        self.calculate_inv_lattice()

    def calc_absolute_coordinates(self, repeat=(1, 1, 1), offset=(0,0,0),edges=False):
        """Returns the cartesion coordinates

        :param repeat: 3 element tuple which determines the number of repeated unit cells
        :param offset: optional offset from (0,0,0)
        :param edges: Determines whether atoms on the faces or edges should be plotted only ones (->False) or multiple times (->True)
        :return: cartesian coordinates of the atoms
        :rtype: Nx4 numpy array
        """
        n_repeat = repeat[0] * repeat[1] * repeat[2]

        abs_coord = np.zeros((self.n_atoms, repeat[0], repeat[1], repeat[2], 4))

        for j1 in range(repeat[0]):
            for j2 in range(repeat[1]):
                for j3 in range(repeat[2]):
                    for i in range(self.n_atoms):
                        abs_coord[i, j1, j2, j3, :3] = np.dot(self.lattice_vectors.T,
                        self.atoms[i, :3].T) + (j3-offset[2]) * self.lattice_vectors[2,:] + (j2-offset[1]) * self.lattice_vectors[1,:] + (j1-offset[0]) * self.lattice_vectors[0, :]
                        abs_coord[i, j1, j2, j3, 3] = self.atoms[i, 3]
                abs_coord_out = abs_coord.reshape((n_repeat * self.n_atoms, 4))
        if edges:
            abs_coord_out = self._add_edges(abs_coord_out,repeat)
        return abs_coord_out

    def find_bonds(self, abs_coords):
        """Searches for bonds between all atoms. This function uses literature values for the covalent bond radii.

        :param abs_coords: cartesion coordinates of the atoms (Nx4) as returned by calc_absolute_coordinates
        :return: Tuple of all bonds between the atoms by their index.
        :rtype: A tuple of two-tuples, e.g. ((3,5),(1,2))

        """
        n_atoms = abs_coords.shape[0]
        abs_coords_pure = abs_coords[:, :3]

        bonds = set()
        cov_radii_array = np.array([cov_radii[int(x)] for x in abs_coords[:, 3]])
        for j1 in range(n_atoms):
            dist = np.linalg.norm(abs_coords_pure[j1 + 1:, :] - abs_coords_pure[j1, :], axis=1)
            mask_array = dist < (cov_radii_array[j1 + 1:] + cov_radii[int(abs_coords[j1, 3])]) * 1.3
            indexes = np.where(mask_array)[0]
            for index in indexes:
                if j1 != index + j1 + 1:
                    bonds.add((j1, index + j1 + 1))

        return tuple(bonds)

    def convert_to_tpiba(self, band_structure_points):
        if type(band_structure_points) in [list, tuple]:
            band_structure_points = np.array(band_structure_points)
        N = band_structure_points.shape[0]
        conv_points = np.zeros((N, 3))
        a = np.linalg.norm(self.lattice_vectors[0, :])
        for i in range(N):
            conv_points[i, :] = np.dot(self.inv_lattice_vectors.T, band_structure_points[i, :]) / (2 * np.pi / a)
        return conv_points

    def calculate_inv_lattice(self):
        volume = np.dot(np.cross(self.lattice_vectors[0, :], self.lattice_vectors[1, :]), self.lattice_vectors[2, :])
        if volume == 0:
            raise ValueError('invalid structre. Two vectors are parallel')
        self.inv_lattice_vectors = np.zeros((3, 3))
        self.inv_lattice_vectors[0, :] = np.cross(self.lattice_vectors[1, :],
                                                  self.lattice_vectors[2, :]) * 2 * np.pi / volume
        self.inv_lattice_vectors[1, :] = np.cross(self.lattice_vectors[2, :],
                                                  self.lattice_vectors[0, :]) * 2 * np.pi / volume
        self.inv_lattice_vectors[2, :] = np.cross(self.lattice_vectors[0, :],
                                                  self.lattice_vectors[1, :]) * 2 * np.pi / volume

    def _add_edges(self,coords,repeat):

        relative_coords = np.zeros(coords.shape)
        for i,coord in enumerate(coords):
            relative_coords[i,:3] = np.dot(np.linalg.inv(self.lattice_vectors.T),coord[:3])
            relative_coords[i,3] = coord[3]

        new_coords = []
        for i,coord in enumerate(relative_coords):
            sym_index = []
            for j,el in enumerate(coord[:3]):
                if abs(el)<1e-12 or abs(el-repeat[j]+1)<1e-12:
                    sym_index.append(j)

            if len(sym_index)>0: # add single vectors
                for j in sym_index:
                    rub = np.zeros((4))
                    rub[:3] = np.dot(self.lattice_vectors.T,coord[:3]) + self.lattice_vectors[j,:]
                    rub[3] = coord[3]
                    new_coords.append(rub)
            if len(sym_index) > 1: # add sum of two vectors
                for pair in itertools.product(sym_index, repeat=2):
                    if pair[0] == pair[1]:
                        continue
                    rub = np.zeros((4))
                    rub[:3] = np.dot(self.lattice_vectors.T,coord[:3]) + self.lattice_vectors[pair[0],:]+self.lattice_vectors[pair[1],:]
                    rub[3] = coord[3]
                    new_coords.append(rub)
            if len(sym_index) > 2: # add sum of three vectors
                rub = np.zeros((4))
                rub[:3] = np.dot(self.lattice_vectors.T, coord[:3]) + self.lattice_vectors[0,:] + self.lattice_vectors[1,:]+ self.lattice_vectors[2,:]
                rub[3] = coord[3]
                new_coords.append(rub)

        new_coords_arr = np.array(new_coords)
        if len(new_coords_arr)>0:
            new_coords_out = np.concatenate((coords,new_coords_arr),axis=0)
        else: new_coords_out = coords

        return new_coords_out

    def density(self,unit='atomic'):
        volume = np.dot(np.cross(self.lattice_vectors[0, :], self.lattice_vectors[1, :]), self.lattice_vectors[2, :])
        species = self.atoms[:,3]

        mass_list = [elemental_masses[int(x)] for x in species]
        mass = sum(mass_list)
        if unit == 'atomic':
            dens = mass/volume
        elif unit == 'g/cm^3':
            dens = mass*unit_mass/(volume*(bohr*1e-8)**3)*1e3
        return dens

    def lattice_information(self):

        lattice = mg.Lattice(self.lattice_vectors)
        atoms = self.atoms
        if len(atoms)==0:
            return {'space group':'','point group':'','crystal system':''}
        structure_mg = mg.Structure(lattice, atoms[:, 3], atoms[:, :3])
        analyzer = SpacegroupAnalyzer(structure_mg)
        res = {'space group':analyzer.get_space_group_symbol(),'point group':analyzer.get_point_group_symbol(),'crystal system':analyzer.get_crystal_system()}
        return res


class BandStructure(object):
    def __init__(self, bands, special_k_points=None, bs_type='electronic'):
        self.bands = bands
        try:
            self.bandgap, self.k_bandgap = self._find_bandgap(bands)
        except Exception:
            self.bandgap, self.k_bandgap = (0, None)
        self.special_k_points = special_k_points
        self.bs_type = bs_type
        self.engine_information = None

    def _find_bandgap(self, bands):
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


class VibrationalStructure(object):
    def __init__(self,frequencies):
        self.frequencies = frequencies
        self.engine_information = None


class EnergyDiagram(object):
    def __init__(self, energies, labels, occupations=None):
        self.energies = energies
        self.labels = labels
        self.occupations = occupations
        self.homo_lumo_gap, self.E_fermi = self._find_homo_lumo_gap(energies)
        self.engine_information = None

    def _find_homo_lumo_gap(self, energies):
        if self.occupations is None:
            index = bisect(energies, 0)
            gap = energies[index] - energies[index - 1]
            E_fermi = energies[index - 1] + gap / 2
        else:
            unoccupied_energies = []
            occupied_energies = []
            for energy, occupation in zip(self.energies, self.occupations):
                if occupation == 0:
                    unoccupied_energies.append(energy)
                else:
                    occupied_energies.append(energy)
            gap = min(unoccupied_energies) - max(occupied_energies)
            E_fermi = max(occupied_energies) + gap / 2

        return gap, E_fermi


class OpticalSpectrum:
    def __init__(self, energy, epsilon2, epsilon1=None):
        self.energy = energy  # energy in eV
        self.engine_information = None
        self.epsilon2_11 = None
        self.epsilon2_22 = None
        self.epsilon2_33 = None
        self.epsilon1_11 = None
        self.epsilon1_22 = None
        self.epsilon1_33 = None

        if type(epsilon2) == list or type(epsilon2) == tuple:
            self.epsilon2_11 = epsilon2[0]
            self.epsilon2_22 = epsilon2[1]
            self.epsilon2_33 = epsilon2[2]
            self.epsilon2 = (self.epsilon2_11 + self.epsilon2_22 + self.epsilon2_33) / 3
        else:
            self.epsilon2 = epsilon2

        if epsilon1 is not None and (type(epsilon1) == list or type(epsilon1) == tuple):
            self.epsilon1_11 = epsilon1[0]
            self.epsilon1_22 = epsilon1[1]
            self.epsilon1_33 = epsilon1[2]
            self.epsilon1 = (self.epsilon1_11 + self.epsilon1_22 + self.epsilon1_33) / 3
        else:
            self.epsilon1 = epsilon1

        self.all_epsilons = [self.epsilon1, self.epsilon1_11, self.epsilon1_22, self.epsilon1_33, self.epsilon2, self.epsilon2_11,
                             self.epsilon1_22, self.epsilon1_33]


class DensityOfStates(object):
    def __init__(self,data,proj_dos=None):
        self.array = np.array(data).astype(np.float)
        self.engine_information = None
        self.proj_dos = proj_dos


class PhononEigenvectors(object):
    def __init__(self,frequencies,modes,k_vectors):
        if frequencies.shape != modes.shape[0:2] or len(modes.shape)!=3 or len(frequencies.shape)!=2:
            raise ValueError('bad shape for input')

        self.frequencies = frequencies
        self.modes = modes
        self.k_vectors = k_vectors

    def calculate_mode(self,n,k,repeat=(1,1,1)):
        n_repeat = repeat[0] * repeat[1] * repeat[2]
        n_atoms = len(self.modes[n,k])//3
        mode = self.modes[n, k]
        mode_conv = mode.reshape((n_atoms,3))
        repeated_mode = np.zeros((n_atoms,repeat[0], repeat[1], repeat[2], 3))

        for j1 in range(repeat[0]):
            for j2 in range(repeat[1]):
                for j3 in range(repeat[2]):
                    for i in range(n_atoms):
                        repeated_mode[i, j1, j2, j3, :] = mode_conv[i,:]

        mode_out = repeated_mode.reshape((n_repeat*n_atoms ,3))
        return mode_out


class StructureParser:
    def __init__(self):
        pass

    def parse_cif_file(self, filename):
        with open(filename, 'r') as f:
            cif_text = f.read()
        cif_text = cif_text.replace('\r', '')
        alpha = float(re.findall(r"_cell_angle_alpha[\s\t]*[-+]?\d*[\.\d+]?", cif_text)[0].split()[1])
        beta = float(re.findall(r"_cell_angle_beta[\s\t]*[-+]?\d*[\.\d+]?", cif_text)[0].split()[1])
        gamma = float(re.findall(r"_cell_angle_gamma[\s\t]*[-+]?\d*[\.\d+]?", cif_text)[0].split()[1])
        a = float(re.findall(r"_cell_length_a[\s\t]*([-+]?\d*\.\d*)", cif_text)[0]) / bohr
        b = float(re.findall(r"_cell_length_b[\s\t]*([-+]?\d*\.\d*)", cif_text)[0]) / bohr
        c = float(re.findall(r"_cell_length_c[\s\t]*([-+]?\d*\.\d*)", cif_text)[0]) / bohr
        params = [a, b, c, alpha, beta, gamma]

        unit_vectors = np.array(calculate_lattice_vectors_from_parameters(params))

        atom_lines, start_of_frac = self.find_atom_lines(cif_text)

        n_atoms = len(atom_lines)
        atom_array = np.zeros((n_atoms, 4))

        for i, atom_line in enumerate(atom_lines):
            atom_l = atom_line.split()
            species = atom_l[0]
            species = self.remove_numbers_from_string(species)
            species = species.title()

            x = atom_l[start_of_frac].split('(')[0]
            y = atom_l[start_of_frac + 1].split('(')[0]
            z = atom_l[start_of_frac + 2].split('(')[0]
            atom_array[i, :] = np.array([x, y, z, p_table_rev[species]])

        atom_array = atom_array[atom_array[:, 3] != 0, :]
        sym_lines = find_lines_between(cif_text, '_symmetry_equiv_pos_as_xyz', 'loop_', strip=True)
        sym_lines = self.remove_cif_attributes(sym_lines)
        n_sym = len(sym_lines)
        sym_atom_array = np.zeros((n_atoms * n_sym, 4))

        if len(sym_lines) > 0:
            counter = 0
            if sym_lines[0][0].isdigit():
                sym_enumeration = True
            else:
                sym_enumeration = False

            for sym_line in sym_lines:
                if sym_enumeration:
                    sym_line = self.remove_counter(sym_line)
                sym = sym_line.replace("'", '').split(',')
                for i in range(n_atoms):
                    pos = atom_array[i, :3]
                    atom_spec = atom_array[i, 3]
                    new_pos = self.perform_sym(pos, sym)
                    sym_atom_array[counter, :3] = new_pos
                    sym_atom_array[counter, 3] = atom_spec
                    counter += 1

            atom_array_finally = remove_duplicates(sym_atom_array)
        else:
            atom_array_finally = atom_array

        atom_array_finally_sorted = atom_array_finally[np.argsort(atom_array_finally[:, 3]), :]

        return CrystalStructure(unit_vectors, atom_array_finally_sorted, scale=a)

    def perform_sym(self, pos, sym):
        x = pos[0]
        y = pos[1]
        z = pos[2]
        namespace = {'x': x, 'y': y, 'z': z}
        exec('xn = ' + sym[0], namespace)
        exec('yn = ' + sym[1], namespace)
        exec('zn = ' + sym[2], namespace)
        xn = namespace["xn"] % 1
        yn = namespace["yn"] % 1
        zn = namespace["zn"] % 1
        return np.array([xn, yn, zn])

    def remove_numbers_from_string(self, x):
        new_string = ''
        for i, el in enumerate(x):
            if not el.isalpha():
                break
            else:
                new_string += el
        return new_string

    def remove_counter(self, x):
        for i, el in enumerate(x):
            if not el.isdigit():
                break
        return x[i:]

    def find_atom_lines(self, text):
        text_lines = text.split('\n')
        text_lines = [x.strip() for x in text_lines]
        for i, line in enumerate(text_lines):
            if line == '_atom_site_fract_x':
                line_of_x = i
                break
        for i in range(1, line_of_x):
            back_line = text_lines[line_of_x - i]
            if back_line[0] != '_':
                start_of_block = line_of_x - i + 1
                break
        number_of_other_lines = line_of_x - start_of_block
        atom_lines = []
        for line in text_lines[line_of_x:]:
            if len(line) == 0 or line == 'loop_':
                break
            if line[0] == '_':
                continue
            else:
                atom_lines.append(line)
        return atom_lines, number_of_other_lines

    def remove_cif_attributes(self, text):
        res = []
        for line in text:
            if line[0] == '_':
                break
            else:
                res.append(line)
        return res


class ComputationalMethods(object):
    def __init__(self, methods):
        self.all_methods = ['periodic', 'non-periodic', 'scf', 'g0w0', 'optical spectrum', 'phonons', 'relax',
                       'bandstructure','dos']
        self.descriptions = {'periodic': 'Periodic structures (crystals)',
                             'non-periodic': 'Non periodic structures (molecules)',
                             'scf': 'Ground state properties', 'g0w0': 'The gw method for many body-corrections',
                             'optical spectrum': 'Calculation of optical spectra',
                             'phonons': 'Phonon properties', 'relax': 'Structure relaxation',
                             'bandstructure': 'Calculation of the band structure','dos':'Calculation of the density of states'}

        if methods is None:
            methods = self.all_methods
        for method in methods:
            if method not in self.all_methods:
                raise ValueError(
                    'methods: ' + str(method) + ' is not known. Available methods are ' + ', '.join(self.all_methods))
        self.methods = methods

    def add(self,method):
        if method not in self.all_methods:
            raise ValueError(
                'methods: ' + str(method) + ' is not known. Available methods are ' + ', '.join(self.all_methods))

        self.methods.append(method)


    def get_description(self, item):
        return self.descriptions[item]

    def __getitem__(self, item):
        return self.methods[item]

    def __setitem__(self, key, value):
        raise TypeError('ComputationalMethods object does not support item assignment')

    def __iter__(self):
        return iter(self.methods)


class GeneralHandler:
    def __init__(self):
        from exciting_handler import Handler
        self.exciting_handler = Handler
        from quantum_espresso_handler import Handler
        self.quantum_espresso_handler = Handler
        from nwchem_handler import Handler
        self.nwchem_handler = Handler
        from abinit_handler import Handler
        self.abinit_handler = Handler
        from ocean_handlers import OceanQe,OceanAbinit,FeffAbinit
        self.ocean_qe_handler = OceanQe
        self.ocean_abi_handler = OceanAbinit
        self.feff_abi_handler = FeffAbinit

        d = {'exciting': self.exciting_handler, 'quantum espresso': self.quantum_espresso_handler,
             'nwchem': self.nwchem_handler, 'abinit': self.abinit_handler,
             "abinit + ocean":self.ocean_abi_handler,'abinit + feff':self.feff_abi_handler}
        self.handlers = OrderedDict(sorted(d.items(), key=lambda t: t[0]))

    def is_handler_available(self, engine_name):
        handler = self.handlers[engine_name]()

        if len(handler.dft_installation_folder) > 0:
            return True
        else:
            return False

    def parse_input_file(self, engine_name, filename):
        handler = self.handlers[engine_name]()
        return handler.parse_input_file(filename)


class KohnShamDensity(object):
    def __init__(self, density):
        self.density = density
        self.engine_information = None


class MolecularDensity(object):
    def __init__(self, density, lattice_vecs, origin):
        self.density = density
        self.grid_vectors = lattice_vecs
        self.origin = origin
        self.engine_information = None


def calculate_lattice_vectors_from_parameters(parameters):
    a, b, c, alpha, beta, gamma = parameters
    alpha = alpha * np.pi / 180
    beta = beta * np.pi / 180
    gamma = gamma * np.pi / 180

    a1 = a * np.array([1, 0, 0])
    a2 = b * np.array([np.cos(gamma), np.sin(gamma), 0])
    x = np.cos(beta)
    y = (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
    z = np.sqrt(1 - x ** 2 - y ** 2)
    a3 = c * np.array([x, y, z])
    return a1, a2, a3


def find_lines_between(text, a, b, strip=False):
    reading = False
    atom_lines = []
    for line in text.split('\n'):
        if strip:
            line = line.strip()
        if line == a:
            reading = True
            continue
        elif reading and line == b:
            break
        if reading and len(line) > 0:
            atom_lines.append(line)
    return atom_lines


def construct_brillouin_vertices(crystal_structure):
    inv_lattice = crystal_structure.inv_lattice_vectors
    l1 = inv_lattice[0,:]
    l2 = inv_lattice[1,:]
    l3 = inv_lattice[2,:]

    origin = 0 * l1
    point_array = np.zeros((3**3, 3))

    counter = 0
    for i in range(-1, 2):
        for j in range(-1, 2):
            for k in range(-1, 2):
                point_array[counter, :] = l1 * i + l2 * j + k * l3
                counter += 1

    all_points = np.append(np.array([origin]), point_array, axis=0)
    voronoi = Voronoi(all_points)
    wigner_points = voronoi.vertices

    wigner_points_cleaned = []

    for w_point in wigner_points:
        dist = np.linalg.norm(point_array - w_point, axis=1)
        if np.all(np.linalg.norm(w_point - origin) <= dist * 1.00001):
            wigner_points_cleaned.append(w_point)

    vertices_array = np.zeros((len(wigner_points_cleaned), 3))
    for i, w_point in enumerate(wigner_points_cleaned):
        vertices_array[i, :] = w_point

    return remove_duplicates(vertices_array)


def construct_convex_hull(w_points):
    hull = ConvexHull(w_points)
    return hull.simplices


def convert_hs_path_to_own(hs_path):
    kpoints = hs_path.kpath['kpoints']
    path = hs_path.kpath['path']

    def convert_path(path, kpoints):
        conv_path = []
        for pos in path[0]:
            if pos == u'\\Gamma':
                pos_r = 'Gamma'
            else:
                pos_r = pos
            new_point = [kpoints[pos], pos_r]
            conv_path.append(new_point)

        return conv_path

    conv_path = convert_path(path, kpoints)
    return conv_path


def calculate_standard_path(structure):
    lattice = mg.Lattice(structure.lattice_vectors)
    atoms = structure.atoms
    structure_mg = mg.Structure(lattice, atoms[:, 3], atoms[:, :3])
    hs_path = HighSymmKpath(structure_mg, symprec=0.1)

    if not hs_path.kpath: # fix for bug in pymatgen that for certain structures no path is returned
        raise Exception('High symmetry path generation failed.')

    conv_path = convert_hs_path_to_own(hs_path)
    return conv_path


def get_emergency_path():
    trash_bs_points = np.array([[0, 0, 0],[0.5, 0, 0],[0, 0.5, 0],[0, 0, 0],[0, 0, 0.5] ])
    trash_bs_labels = ['GAMMA','X','Y','GAMMA','Z']
    path = list(zip(trash_bs_points, trash_bs_labels))
    return path


def calculate_path_length(structure, k_path):
    points = []
    pos = 0
    p_length = 0
    last_point = None
    for point, label in k_path:
        if last_point is None:
            points.append([pos, label])
            last_point = np.dot(structure.inv_lattice_vectors, point)
        else:
            abs_point = np.dot(structure.inv_lattice_vectors, point)
            p_length = np.linalg.norm(abs_point - last_point)
            points.append([pos + p_length, label])
            last_point = abs_point
        pos += p_length
    return points


def bonds_to_path(bonds_in):
    bonds = copy.deepcopy(bonds_in)
    paths = []
    path = []

    while len(bonds) > 0:

        if len(path) == 0:
            path = list(bonds[0])
            del bonds[0]
        else:
            last_el = path[-1]
            for i, bond in enumerate(bonds):
                if last_el == bond[0]:
                    path.append(bond[1])
                    del bonds[i]
                    break
                elif last_el == bond[1]:
                    path.append(bond[0])
                    del bonds[i]
                    break
                elif i == len(bonds) - 1:
                    paths.append(path)
                    path = []
    if len(path) > 0:
        paths.append(path)
    return paths


def calculate_mass_matrix(structure,repeat=(1,1,1),edges=False,phonon_conv=False):
    coords = structure.calc_absolute_coordinates(repeat=repeat,edges=edges)
    species = coords[:,3]
    n_repeat = repeat[0]*repeat[1]*repeat[2]

    mass_matrix = np.zeros((coords.shape[0],coords.shape[0]))
    if phonon_conv:
        mass_list = np.array([[np.sqrt(elemental_masses[int(x)])**-1] for x in species ])
    else:
        mass_list = np.array([[elemental_masses[int(x)]] for x in species ])
    mass_list_coll = [item for sub_list in mass_list for item in sub_list]
    np.fill_diagonal(mass_matrix,mass_list_coll)
    return mass_matrix


def convert_pymatgen_structure(structure):
    cart_coords = structure.cart_coords
    atoms_cart = np.zeros((cart_coords.shape[0],4))
    atoms_cart[:,:3] = structure.cart_coords/bohr
    for i,specie in enumerate(structure.species):
        atoms_cart[i,3] = specie.Z

    crystal_structure = CrystalStructure(structure.lattice.matrix/bohr, atoms_cart,relative_coords=False)
    return crystal_structure


def query_materials_database(structure_formula):
    with MPRester(API_KEY) as m:
        data = m.get_data(structure_formula)
    return data


def get_materials_structure_from_id(id):
    with MPRester(API_KEY) as m:
        structure = m.get_structure_by_material_id(id)
    return convert_pymatgen_structure(structure)


def replace_greek(x):
    if x == '\\Gamma':
        res = 'Gamma'
    else:
        res = x
    return res


def find_center_of_faces(coordinates, triangles):
    faces = construct_faces(coordinates, triangles)
    centers = np.zeros((len(faces), 3))

    for i, vertex in enumerate(faces):
        centers[i, :] = np.mean(coordinates[list(vertex), :], axis=0)
    return centers


def is_in_plane(base_triangle, test_triangle, coordinates):
    normal = construct_normal_from_triangle(coordinates, base_triangle)
    normal_alt = construct_normal_from_triangle(coordinates, test_triangle)

    if base_triangle[0] == test_triangle[0]:
        test_vec = coordinates[base_triangle[0], :] - coordinates[test_triangle[1], :]
    else:
        test_vec = coordinates[base_triangle[0], :] - coordinates[test_triangle[0], :]

    test_vec = test_vec/np.linalg.norm(test_vec)
    return (abs(np.dot(normal, normal_alt)) > 0.999) and (abs(np.dot(normal, test_vec)) < 0.001)


def construct_normal_from_triangle(coordinates, triangle):
    v1 = coordinates[triangle[1], :] - coordinates[triangle[0], :]
    v2 = coordinates[triangle[2], :] - coordinates[triangle[0], :]
    normal = np.cross(v1, v2)
    return normal/np.linalg.norm(normal)


def is_list_in_list_of_lists(small_list, large_list):
    for el in large_list:
        if sorted(small_list) == sorted(el):
            return True
    return False


def construct_faces(points, triangles):
    used_triangles = []
    faces = []
    for base_triangle in triangles:
        if is_list_in_list_of_lists(base_triangle, used_triangles):
            continue
        used_triangles.append(base_triangle)
        face = set(base_triangle)
        for test_triangle in triangles:
            if is_list_in_list_of_lists(test_triangle, used_triangles):
                continue
            if is_in_plane(base_triangle, test_triangle, points):
                face.update(test_triangle)
                used_triangles.append(test_triangle)
        faces.append(list(face))

    return faces


def get_materials_dos_from_id(id):
    with MPRester(API_KEY) as m:
        dos = m.get_dos_by_material_id(id)

    data = np.zeros((len(dos.energies),2))
    data[:,0] = dos.energies - dos.get_gap()/2
    for key,value in dos.densities.items():
        data[:,1] = data[:,1] + value
    return DensityOfStates(data)


def get_materials_band_structure_from_id(id):
    with MPRester(API_KEY) as m:
        band_structure = m.get_bandstructure_by_material_id(id)
        structure = m.get_structure_by_material_id(id)
    conv_structure = convert_pymatgen_structure(structure)
    bands = []
    for key,band in band_structure.bands.items():
        bands.append(band)
    conv_bands = []

    k = []
    for kpoint in band_structure.kpoints:
        k.append([kpoint.frac_coords,''])
    k_scalar = np.array(calculate_path_length(conv_structure,k))[:,0]

    for band in bands:
        for i in range(band.shape[0]):
            array = np.zeros((len(band[i,:]),2))
            array[:,0] = k_scalar
            array[:,1] = band[i,:] - band_structure.efermi-band_structure.get_band_gap()['energy']/2
            conv_bands.append(array)

    sorted_bands = conv_bands
    # sorted_bands = sorted(conv_bands,key=np.mean)
    # high_sym_path = HighSymmKpath(structure)
    # path = convert_hs_path_to_own(high_sym_path)
    path = []
    last_label = None
    for branch in band_structure.branches:
        s_index = branch['start_index']
        e_index = branch['end_index']

        label = replace_greek(branch['name'].split('-')[0])

        if last_label is not None and last_label!= label:
            path.append([last_coords, last_label])
        last_label = replace_greek(branch['name'].split('-')[1])
        last_coords = band_structure.kpoints[e_index].frac_coords

        path.append([band_structure.kpoints[s_index].frac_coords,label])


    e_index = branch['end_index']
    label = branch['name'].split('-')[1]
    if label == u'\\Gamma':
        label = 'Gamma'
    path.append([band_structure.kpoints[e_index].frac_coords, label])

    path_conv = calculate_path_length(conv_structure,path)
    return BandStructure(sorted_bands,special_k_points=path_conv)


if __name__ == "__main__":
    data = query_materials_database('LiBH4')
    id = data[0]['material_id']

    id = 'mp-66'
    structure = get_materials_structure_from_id(id)
    # bs = get_materials_band_structure_from_id(id)
    dos = get_materials_dos_from_id(id)
    coords = structure.calc_absolute_coordinates()



    # import mayavi.mlab as mlab
    #
    # mlab.points3d(coords[:, 0], coords[:, 1], coords[:, 2], color=(0.7, 0.7, 0.7), scale_factor=1, )
    # mlab.show()


    # atoms = np.array([[0, 0, 0, 6], [0.333333333333333, 0.3333333333333333333, 0.0, 10]])
    # unit_cell = 4.650000* np.array([[0.5, 0.866025, 0], [-0.5, 0.866025, 0.0], [0, 0.0, 6.0]])
    # # unit_cell = 6.719 * np.array([[0.5, 0.5, -0.5], [0.5, -0.5, 0.5], [-0.5, 0.5, 0.5]])
    # #
    #
    #
    # # phi = 60/180*np.pi
    # # unit_cell = 6.719 * np.array([ [1,0,0], [np.cos(phi),np.sin(phi),0], [0, 0, 3]])
    #
    # p = StructureParser()
    # data = p.parse_cif_file('/home/jannick/OpenDFT_projects/libh4_theo/theo.cif')
    #
    # crystal_structure = CrystalStructure(unit_cell, atoms)
    # mass_matrix = calculate_mass_matrix(crystal_structure)
    # coords = crystal_structure.calc_absolute_coordinates(repeat=[2, 1, 1],edges=True)
    #
    # density = crystal_structure.density(unit='g/cm^3')
    # print(density)

    # bonds = crystal_structure.find_bonds(coords)
    # print(crystal_structure.find_bonds(coords))
    #
    # path = calculate_standard_path(crystal_structure)
    #
    # from visualization import BrillouinVisualization
    #
    # vis = BrillouinVisualization(None)
    # vis.set_crystal_structure(crystal_structure)
    # vis.set_path(path )
    # vis.configure_traits()

    # import mayavi.mlab as mlab
    #
    # ts = time.time()
    # w_points = construct_brillouin_vertices(crystal_structure)
    # brillouin_edges = construct_convex_hull(w_points)
    #
    #
    # mlab.points3d(w_points[:,0],w_points[:,1],w_points[:,2],color=(0.7, 0.7, 0.7),scale_factor=.1,)
    #
    # mlab.triangular_mesh(w_points[:,0],w_points[:,1],w_points[:,2],brillouin_edges,opacity=0.7,color=(0.5,0.5,0.5))
    #
    # # for i,connection in enumerate(brillouin_edges):
    # #     for con in connection:
    # #         bond = [i,con]
    # #         mlab.plot3d(w_points[bond, 0], w_points[bond, 1], w_points[bond, 2])
    #
    # mlab.show()


    # parser = StructureParser()
    # out = parser.parse_cif_file('/home/jannick/OpenDFT_projects/LiBH4/1504402.cif')
    # general_handler = GeneralHandler()
    # print(general_handler.available_handlers)


