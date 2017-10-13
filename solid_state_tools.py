from __future__ import division
import numpy as np
import re
import periodictable as pt
from bisect import bisect
import time
from scipy.spatial import ConvexHull,Voronoi

bohr = 0.52917721067

cov_radii = np.loadtxt('./data/cov_radii.dat')/bohr

p_table = {i: el.__repr__() for i, el in enumerate(pt.elements)}
p_table_rev = {el.__repr__(): i for i, el in enumerate(pt.elements)}


def remove_duplicates_old(data, treshold=0.01):
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
    def __init__(self, atoms,scale=1.0):
        self.atoms = np.array(atoms,dtype=np.float) # np array with [x,y,z,type] type is number in periodic system
        self.atoms[:,:3] = self.atoms[:,:3]*scale
        self.n_atoms = atoms.shape[0]
        self.scale = scale  # This is just bonus info. Do not use this here. Only for editing

    def calc_absolute_coordinates(self,repeat=[1,1,1]):
        return self.atoms

    def find_bonds(self,abs_coords):
        n_atoms = abs_coords.shape[0]
        abs_coords_pure = abs_coords[:,:3]

        dist_mat = np.zeros((n_atoms,n_atoms))
        bonds = []
        for j1 in range(n_atoms):
            for j2 in range(j1+1,n_atoms):

                dist = np.linalg.norm(abs_coords_pure[j1,:]-abs_coords_pure[j2,:] )
                if dist < (cov_radii[int(abs_coords[j1,3])]+cov_radii[int(abs_coords[j2,3])])*1.3:
                    dist_mat[j1, j2] = dist
                    bonds.append([j1,j2]);
        return bonds

class CrystalStructure(object):

    def __init__(self,lattice_vectors,atoms,relative_coords=True,scale=1.0):
        self.lattice_vectors = np.array(lattice_vectors,dtype=np.float) # tuple of np.arrays
        volume = np.dot(np.cross(self.lattice_vectors[0,:],self.lattice_vectors[1,:]),self.lattice_vectors[2,:])
        self.inv_lattice_vectors = np.zeros((3,3))
        self.inv_lattice_vectors[0,:] = np.cross(self.lattice_vectors[1,:],self.lattice_vectors[2,:])*2*np.pi/volume
        self.inv_lattice_vectors[1,:] = np.cross(self.lattice_vectors[2,:],self.lattice_vectors[0,:])*2*np.pi/volume
        self.inv_lattice_vectors[2,:] = np.cross(self.lattice_vectors[0,:],self.lattice_vectors[1,:])*2*np.pi/volume

        self.atoms = np.array(atoms,dtype=np.float) # np array with [x,y,z,type] type is number in periodic system
        self.n_atoms = atoms.shape[0]
        self.scale = scale  # This is just bonus info. Do not use this here. Only for editing

        if relative_coords:
            self.atoms[:,:3] = np.mod(self.atoms[:,:3],1)
        else:
            inv_lattice_vecs = np.linalg.inv(self.lattice_vectors.T)
            for i in range(self.n_atoms):
                self.atoms[i,:3] = np.dot(inv_lattice_vecs,self.atoms[i,:3].T)

    def calc_absolute_coordinates(self,repeat=[1,1,1]):
        n_repeat = repeat[0]*repeat[1]*repeat[2]

        abs_coord = np.zeros((self.n_atoms,repeat[0],repeat[1],repeat[2],4))

        for j1 in range(repeat[0]):
            for j2 in range(repeat[1]):
                for j3 in range(repeat[2]):
                    for i in range(self.n_atoms):
                        abs_coord[i,j1,j2,j3,:3] = np.dot(self.lattice_vectors.T,self.atoms[i,:3].T) + j3*self.lattice_vectors[2,:]+j2*self.lattice_vectors[1,:]+j1*self.lattice_vectors[0,:]
                        abs_coord[i,j1,j2,j3,3] = self.atoms[i,3]
                abs_coord_out = abs_coord.reshape((n_repeat*self.n_atoms,4))
        return abs_coord_out

    def find_bonds(self,abs_coords):
        n_atoms = abs_coords.shape[0]
        abs_coords_pure = abs_coords[:,:3]

        dist_mat = np.zeros((n_atoms,n_atoms))
        bonds = []
        for j1 in range(n_atoms):
            for j2 in range(j1+1,n_atoms):

                dist = np.linalg.norm(abs_coords_pure[j1,:]-abs_coords_pure[j2,:] )
                if dist < (cov_radii[int(abs_coords[j1,3])]+cov_radii[int(abs_coords[j2,3])])*1.3:
                    dist_mat[j1, j2] = dist
                    bonds.append([j1,j2]);
        return bonds

    def convert_to_tpiba(self,band_structure_points):
        if type(band_structure_points) in [list,tuple]:
            band_structure_points = np.array(band_structure_points)
        N = band_structure_points.shape[0]
        conv_points = np.zeros((N,3))
        a = np.linalg.norm(self.lattice_vectors[0, :])
        for i in range(N):
            conv_points[i,:] = np.dot(self.inv_lattice_vectors.T,band_structure_points[i,:])/(2*np.pi/a)
        return conv_points

class BandStructure(object):
    def __init__(self,bands,special_k_points=None,bs_type='electronic'):
        self.bands = bands
        self.bandgap,self.k_bandgap = self._find_bandgap(bands)
        self.special_k_points = special_k_points
        self.bs_type = bs_type

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

class EnergyDiagram(object):
    def __init__(self,energies,labels,occupations=None):
        self.energies = energies
        self.labels = labels
        self.occupations = occupations
        self.homo_lumo_gap,self.E_fermi = self._find_homo_lumo_gap(energies)


    def _find_homo_lumo_gap(self,energies):
        if self.occupations is None:
            index = bisect(energies,0)
            gap = energies[index] - energies[index-1]
            E_fermi = energies[index-1] + gap/2
        else:
            unoccupied_energies = []
            occupied_energies = []
            for energy,occupation in zip(self.energies,self.occupations):
                if occupation==0:
                    unoccupied_energies.append(energy)
                else:
                    occupied_energies.append(energy)
            gap = min(unoccupied_energies) - max(occupied_energies)
            E_fermi = max(occupied_energies) + gap/2

        return gap,E_fermi


class OpticalSpectrum:
    def __init__(self,energy,epsilon2,epsilon1=None):
        self.energy = energy # energy in eV

        if type(epsilon2) == list or type(epsilon2)== tuple:
            self.epsilon2_11 = epsilon2[0]
            self.epsilon2_22 = epsilon2[1]
            self.epsilon2_33 = epsilon2[2]
            self.epsilon2 = (self.epsilon2_11+self.epsilon2_22+self.epsilon2_33)/3
        else:
            self.epsilon2 = epsilon2

        if epsilon1 is not None and (type(epsilon1) == list or type(epsilon1)== tuple):
            self.epsilon1_11 = epsilon1[0]
            self.epsilon1_22 = epsilon1[1]
            self.epsilon1_33 = epsilon1[2]
            self.epsilon1 = (self.epsilon1_11 + self.epsilon1_22 + self.epsilon1_33) / 3
        else:
            self.epsilon1 = epsilon1

class StructureParser:
    def __init__(self):
        pass

    def parse_cif_file(self,filename):
        with open(filename, 'r') as f:
            cif_text = f.read()
        cif_text = cif_text.replace('\r','')
        alpha = float(re.findall(r"_cell_angle_alpha[\s\t]*[-+]?\d*[\.\d+]?", cif_text)[0].split()[1])
        beta = float(re.findall(r"_cell_angle_beta[\s\t]*[-+]?\d*[\.\d+]?", cif_text)[0].split()[1])
        gamma = float(re.findall(r"_cell_angle_gamma[\s\t]*[-+]?\d*[\.\d+]?", cif_text)[0].split()[1])
        a = float(re.findall(r"_cell_length_a[\s\t]*[-+]?\d*[\.\d+]?", cif_text)[0].split()[1])/bohr
        b = float(re.findall(r"_cell_length_b[\s\t]*[-+]?\d*[\.\d+]?", cif_text)[0].split()[1])/bohr
        c = float(re.findall(r"_cell_length_c[\s\t]*[-+]?\d*[\.\d+]?", cif_text)[0].split()[1])/bohr
        params = [a,b,c,alpha,beta,gamma]

        unit_vectors = np.array(calculate_lattice_vectors_from_parameters(params))


        atom_lines,start_of_frac = self.find_atom_lines(cif_text)

        n_atoms = len(atom_lines)
        atom_array = np.zeros((n_atoms,4))

        for i,atom_line in enumerate(atom_lines):
            atom_l = atom_line.split()
            species = atom_l[0]
            species = self.remove_numbers_from_string(species)
            species = species.title()

            x = atom_l[start_of_frac].split('(')[0]
            y = atom_l[start_of_frac+1].split('(')[0]
            z = atom_l[start_of_frac+2].split('(')[0]
            atom_array[i,:] = np.array([x,y,z,p_table_rev[species]])


        atom_array = atom_array[atom_array[:,3]!=0,:]
        sym_lines = find_lines_between(cif_text,'_symmetry_equiv_pos_as_xyz','loop_',strip=True)
        sym_lines = self.remove_cif_attributes(sym_lines)
        n_sym = len(sym_lines)
        sym_atom_array = np.zeros((n_atoms*n_sym,4))

        counter = 0
        if sym_lines[0][0].isdigit():
            sym_enumeration = True
        else:
            sym_enumeration = False

        for sym_line in sym_lines:
            if sym_enumeration:
                sym_line = self.remove_counter(sym_line)
            sym = sym_line.replace("'",'').split(',')
            for i in range(n_atoms):
                pos = atom_array[i,:3]
                atom_spec = atom_array[i,3]
                new_pos = self.perform_sym(pos,sym)
                sym_atom_array[counter,:3] = new_pos
                sym_atom_array[counter,3] = atom_spec
                counter +=1

        atom_array_finally = remove_duplicates_old(sym_atom_array)
        atom_array_finally_sorted = atom_array_finally[np.argsort(atom_array_finally[:,3]),:]
        return CrystalStructure(unit_vectors,atom_array_finally_sorted,scale=a)


    def perform_sym(self,pos,sym):
        x = pos[0]
        y = pos[1]
        z = pos[2]
        namespace = {'x':x,'y':y,'z':z}
        exec('xn = ' + sym[0],namespace)
        exec('yn = ' + sym[1],namespace)
        exec('zn = ' + sym[2],namespace)
        xn = namespace["xn"] % 1
        yn = namespace["yn"] % 1
        zn = namespace["zn"] % 1
        return np.array([xn,yn,zn])

    def remove_numbers_from_string(self,x):
        new_string = ''
        for i,el in enumerate(x):
            if not el.isalpha():
                break
            else:
                new_string += el
        return new_string

    def remove_counter(self,x):
        for i,el in enumerate(x):
            if not el.isdigit():
                break
        return x[i:]

    def find_atom_lines(self,text):
        text_lines = text.split('\n')
        text_lines = [x.strip() for x in text_lines]
        for i,line in enumerate(text_lines):
            if line == '_atom_site_fract_x':
                line_of_x = i
                break
        for i in range(1,line_of_x):
            back_line = text_lines[line_of_x-i]
            if back_line[0] != '_':
                start_of_block = line_of_x-i+1
                break
        number_of_other_lines = line_of_x - start_of_block
        atom_lines = []
        for line in text_lines[line_of_x:]:
            if len(line)==0 or line == 'loop_':
                break
            if line[0] == '_':
                continue
            else:
                atom_lines.append(line)
        return atom_lines,number_of_other_lines

    def remove_cif_attributes(self,text):
        res = []
        for line in text:
            if line[0] == '_':
                break
            else:
                res.append(line)
        return res

class ComputationalMethods(object):
    def __init__(self,methods):
        all_methods = ['periodic','non-periodic','scf','gw','optical spectrum','phonons','relax']

        if methods is None:
            methods = all_methods
        for method in methods:
            if method not in all_methods:
                raise ValueError('methods: ' +str(method)+' is not known. Available methods are '+', '.join(all_methods))
        self.methods = methods

    def __getitem__(self, item):
        return self.methods[item]

    def __setitem__(self, key, value):
        raise TypeError('ComputationalMethods object does not support item assignment')

    def __iter__(self):
        return iter(self.methods)

class GeneralHandler():
    def __init__(self):
        from exciting_handler import Handler
        self.exciting_handler = Handler
        from quantum_espresso_handler import Handler
        self.quantum_espresso_handler = Handler
        from nwchem_handler import Handler
        self.nwchem_handler = Handler

        self.handlers = {'exciting':self.exciting_handler,'quantum espresso': self.quantum_espresso_handler,'nwchem': self.nwchem_handler}


    def is_handler_available(self,engine_name):
        handler = self.handlers[engine_name]()

        if len(handler.dft_installation_folder) > 0:
            return True
        else:
            return False


    def parse_input_file(self,engine_name,filename):
        handler = self.handlers[engine_name]()
        return handler.parse_input_file(filename)


class KohnShamDensity(object):
    def __init__(self,density):
        self.density = density

class MolecularDensity(object):
    def __init__(self,density,lattice_vecs,origin):
        self.density = density
        self.grid_vectors = lattice_vecs
        self.origin = origin

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

def find_lines_between(text,a,b,strip=False):
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
        if reading and len(line)>0:
            atom_lines.append(line)
    return atom_lines

def construct_brillouin_vertices(crystal_structure):
    inv_lattice = crystal_structure.inv_lattice_vectors
    l1 = inv_lattice[:,0]
    l2 = inv_lattice[:,1]
    l3 = inv_lattice[:,2]

    origin = 0 * l1
    point_array = np.zeros((27,3))

    counter= 0
    for i in range(-1, 2):
        for j in range(-1, 2):
            for k in range(-1, 2):
                point_array[counter,:] = l1 * i + l2 * j + k * l3
                counter += 1

    all_points = np.append(np.array([origin]),point_array,axis=0)
    voronoi = Voronoi(all_points)
    wigner_points = voronoi.vertices

    ### Own voronoi implementation (might still be usefull one day...(okay i am just keeping it because it was hard to do))
    # N_points = point_array.shape[0]
    # wigner_points = []
    # x1, y1, z1 = origin
    #
    # for i in range(N_points):
    #     for j in range(i+1,N_points):
    #         for k in range(j+1,N_points):
    #
    #             if i == j or i == k or j == k:
    #                 continue
    #
    #             x2, y2, z2 = point_array[i,:]
    #             x3, y3, z3 = point_array[j,:]
    #             x4, y4, z4 = point_array[k,:]
    #
    #             A = np.array([[2 * x1 - x2 - x3 - x4, 2 * y1 - y2 - y3 - y4, 2 * z1 - z2 - z3 - z4],
    #                           [2 * x2 - x1 - x3 - x4, 2 * y2 - y1 - y3 - y4, 2 * z2 - z1 - z3 - z4],
    #                           [2 * x3 - x1 - x2 - x4, 2 * y3 - y1 - y2 - y4, 2 * z3 - z1 - z2 - z4]
    #                           ])
    #             r = np.linalg.matrix_rank(A)
    #             if r < 3:
    #                 continue
    #             B = np.array([[x1 ** 2 - 0.5 * (x2 ** 2 + x3 ** 2 + x4 ** 2) + y1 ** 2 - 0.5 * (
    #             y2 ** 2 + y3 ** 2 + y4 ** 2) + z1 ** 2 - 0.5 * (z2 ** 2 + z3 ** 2 + z4 ** 2)],
    #                           [x2 ** 2 - 0.5 * (x1 ** 2 + x3 ** 2 + x4 ** 2) + y2 ** 2 - 0.5 * (
    #                           y1 ** 2 + y3 ** 2 + y4 ** 2) + z2 ** 2 - 0.5 * (z1 ** 2 + z3 ** 2 + z4 ** 2)],
    #                           [x3 ** 2 - 0.5 * (x1 ** 2 + x2 ** 2 + x4 ** 2) + y3 ** 2 - 0.5 * (
    #                           y1 ** 2 + y2 ** 2 + y4 ** 2) + z3 ** 2 - 0.5 * (z1 ** 2 + z2 ** 2 + z4 ** 2)]])
    #
    #             xout = np.dot(np.linalg.inv(A), B).T
    #             xout = np.array([xout[0, 0], xout[0, 1], xout[0, 2]])
    #             wigner_points.append(xout)

    wigner_points_cleaned = []

    for w_point in wigner_points:
        dist = np.linalg.norm(point_array-w_point,axis=1)
        if np.all(np.linalg.norm(w_point - origin) <= dist * 1.01):
            wigner_points_cleaned.append(w_point)

    vertices_array = np.zeros((len(wigner_points_cleaned),3))
    for i,w_point in enumerate(wigner_points_cleaned):
        vertices_array[i,:] = w_point

    return remove_duplicates_old(vertices_array)

def construct_convex_hull(w_points):
    hull = ConvexHull(w_points)
    return hull.simplices

    # bonds = []
    # for simplex in hull.simplices:
    #     simplex = np.append(simplex, simplex[0])
    #     for i in range(len(simplex)-1):
    #         bonds.append([simplex[i],simplex[i+1]])
    #
    # connections = []
    # for i in range(w_points.shape[0]):
    #     con_temp = set()
    #     for bond in bonds:
    #         if i==bond[0]:
    #             con_temp.add(bond[1])
    #         elif i==bond[1]:
    #             con_temp.add(bond[0])
    #     connections.append(con_temp)
    #
    # shortest_connections = []
    # for i,connection in enumerate(connections):
    #     connection = list(connection)
    #     dist = np.linalg.norm(w_points[connection,:]-w_points[i,:] ,axis=1)
    #     in_sort = np.argsort(dist)[:3]
    #     shortest_connections.append(np.array(connection)[in_sort])
    # shortest_connections = connections

    # return shortest_connections

if __name__ == "__main__":
    atoms = np.array([[0, 0, 0, 6], [0.25, 0.25, 0.25, 6]])
    unit_cell = 6.719 * np.array([[0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    # unit_cell = 6.719 * np.array([[0.5, 0.5, -0.5], [0.5, -0.5, 0.5], [-0.5, 0.5, 0.5]])
    # unit_cell = 6.719 * np.eye(3,3);    unit_cell[2,1] = 0.5

    # phi = 60/180*np.pi
    # unit_cell = 6.719 * np.array([ [1,0,0], [np.cos(phi),np.sin(phi),0], [0, 0, 3]])


    crystal_structure = CrystalStructure(unit_cell, atoms)
    coords = crystal_structure.calc_absolute_coordinates(repeat=[1,1,2])
    print(crystal_structure.find_bonds(coords))

    import mayavi.mlab as mlab

    ts = time.time()
    w_points = construct_brillouin_vertices(crystal_structure)
    brillouin_edges = construct_convex_hull(w_points)


    mlab.points3d(w_points[:,0],w_points[:,1],w_points[:,2],color=(0.7, 0.7, 0.7),scale_factor=.1,)

    mlab.triangular_mesh(w_points[:,0],w_points[:,1],w_points[:,2],brillouin_edges,opacity=0.7,color=(0.5,0.5,0.5))

    # for i,connection in enumerate(brillouin_edges):
    #     for con in connection:
    #         bond = [i,con]
    #         mlab.plot3d(w_points[bond, 0], w_points[bond, 1], w_points[bond, 2])

    mlab.show()


        # parser = StructureParser()
    # out = parser.parse_cif_file('/home/jannick/OpenDFT_projects/LiBH4/1504402.cif')
    # general_handler = GeneralHandler()
    # print(general_handler.available_handlers)