from pymatgen.symmetry.bandstructure import HighSymmKpath
import pymatgen as mg
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import spglib as spg
import pickle

test = mg.Structure.from_file('low_temp.cif')
finder = SpacegroupAnalyzer(test, 0.001)
hs_path = HighSymmKpath(test,symprec=0.001) # unabgefangener error mit symprec 0.01 funktioniert mit symprec 0.001


with open('save.pkl', 'rb') as handle:
    b = pickle.load(handle)
    structure = b.pop('crystal structure',None)

lattice = mg.Lattice(structure.lattice_vectors)
atoms = structure.atoms
structure_mg = mg.Structure(lattice, atoms[:, 3], atoms[:, :3])

hs_path = HighSymmKpath(structure_mg,symprec=0.001) # symprec 0.001 ergibt eine kubische Struktur

kpoints = hs_path.kpath['kpoints']
path = hs_path.kpath['path']





