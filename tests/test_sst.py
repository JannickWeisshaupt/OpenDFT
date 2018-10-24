import sys, os

myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, myPath + '/../')

os.environ['ETS_TOOLKIT'] = 'qt4'

import imp

try:
    imp.find_module('PySide')  # test if PySide if available
except ImportError:
    os.environ['QT_API'] = 'pyqt'  # signal to pyface that PyQt4 should be used


import src.solid_state_tools as sst
import numpy as np

unit_cell = 6.719 * np.array([[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]])
atoms = np.array([[0.0, 0.0, 0.0, 6], [0.25, 0.25, 0.25, 6]])
diamond = sst.CrystalStructure(unit_cell, atoms)


h2o_atoms = np.array([[0.12892522000,0.12893350000,0 , 8], [1.94137744000, -0.07030354000, 0, 1],[-0.07030266000, 1.94137003000, 0, 1]])
h2o = sst.MolecularStructure(h2o_atoms)

def compare_path(p,o):
    is_equal = True
    for el1,el2 in zip(p,o):
        dist = np.linalg.norm(el1[0]-el2[0])
        if dist > 1e-10 or el1[1] != el2[1]:
            is_equal = False
            break
    return is_equal


def test_molecular_structure():
    all_tests = []

    info = h2o.symmetry_information()
    bonds = h2o.find_bonds(h2o.calc_absolute_coordinates())

    info_test = str(info['point group']) == 'C2v'
    all_tests.append(info_test)

    bond_test = set(bonds) == set(((0, 1), (0, 2)))
    all_tests.append(bond_test)

    assert all(all_tests)


def test_crystal_structure():
    dens = diamond.density()
    all_tests = []
    coords = []
    bonds = []
    repeat = [[1,1,1],[2,2,2],[5,2,1],[2,3,5]]
    for rep in repeat:
        coord = diamond.calc_absolute_coordinates(rep)
        coords.append(coord)
        bonds.append(diamond.find_bonds(coord))

    correct_number_of_bonds = [1,20,23,89]
    number_of_bonds = [len(x) for x in bonds]

    bond_test = number_of_bonds == correct_number_of_bonds
    all_tests.append(bond_test)

    dens_test = np.abs(dens - 0.31677024692785494) < 1e-12
    all_tests.append(dens_test)

    info = diamond.lattice_information()
    diamond_info = {'crystal system': 'cubic', 'point group': 'm-3m', 'space group': 'Fd-3m'}

    info_test = info == diamond_info
    all_tests.append(info_test)

    assert all(all_tests)


def test_standard_path():
    structure = diamond
    path = sst.calculate_standard_path(structure)
    stored_path = [[np.array([0., 0., 0.]), 'Gamma'], [np.array([0.5, 0. , 0.5]), 'X'], [np.array([0.5 , 0.25, 0.75]), 'W'], [np.array([0.375, 0.375, 0.75 ]), 'K'], [np.array([0., 0., 0.]), 'Gamma'], [np.array([0.5, 0.5, 0.5]), 'L'], [np.array([0.625, 0.25 , 0.625]), 'U'], [np.array([0.5 , 0.25, 0.75]), 'W'], [np.array([0.5, 0.5, 0.5]), 'L'], [np.array([0.375, 0.375, 0.75 ]), 'K']]
    assert compare_path(stored_path,path)


def test_brillouin_construction():
    w_points = sst.construct_brillouin_vertices(diamond)
    brillouin_edges = sst.construct_convex_hull(w_points)
    assert (len(w_points) == 24) and (len(brillouin_edges) == 44)

if __name__ == "__main__":
    test_brillouin_construction()
    test_crystal_structure()
    test_molecular_structure()