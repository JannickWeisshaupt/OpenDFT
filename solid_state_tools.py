import numpy as np

from visualization import cov_radii

class CrystalStructure:

    def __init__(self,lattice_vectors,atoms):
        self.lattice_vectors = np.array(lattice_vectors,dtype=np.float) # tuple of np.arrays
        self.atoms = np.array(atoms,dtype=np.float) # np array with [x,y,z,type] type is number in periodic system
        self.atoms[:,:3] = np.mod(self.atoms[:,:3],1)
        self.n_atoms = atoms.shape[0]

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

class BandStructure:
    def __init__(self,bands,bandgap,k_bandgap,special_k_points=None):
        self.bands = bands
        self.bandgap = bandgap
        self.special_k_points = special_k_points
        self.k_bandgap = k_bandgap

if __name__ == "__main__":

    atoms = np.array([[0,0,0,6],[0.25,0.25,0.25,6]])
    unit_cell = 6.719*np.array([[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5]])

    crystal_structure = CrystalStructure(unit_cell,atoms)
    coords = crystal_structure.calc_absolute_coordinates(repeat=[1,1,2])
    print(crystal_structure.find_bonds(coords))