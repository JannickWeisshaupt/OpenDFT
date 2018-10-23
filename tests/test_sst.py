import sys
sys.path.insert(0,'./src/')
import src.solid_state_tools as sst
import numpy as np

unit_cell = 6.719 * np.array([[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]])
atoms = np.array([[0.0, 0.0, 0.0, 6], [0.25, 0.25, 0.25, 6]])
diamond = sst.CrystalStructure(unit_cell, atoms)


def compare_path(p,o):
    is_equal = True
    for el1,el2 in zip(p,o):
        dist = np.linalg.norm(el1[0]-el2[0])
        if dist > 1e-10 or el1[1] != el2[1]:
            is_equal = False
            break
    return is_equal


def test_standard_path():
    structure = diamond
    path = sst.calculate_standard_path(structure)
    stored_path = [[np.array([0., 0., 0.]), 'Gamma'], [np.array([0.5, 0. , 0.5]), 'X'], [np.array([0.5 , 0.25, 0.75]), 'W'], [np.array([0.375, 0.375, 0.75 ]), 'K'], [np.array([0., 0., 0.]), 'Gamma'], [np.array([0.5, 0.5, 0.5]), 'L'], [np.array([0.625, 0.25 , 0.625]), 'U'], [np.array([0.5 , 0.25, 0.75]), 'W'], [np.array([0.5, 0.5, 0.5]), 'L'], [np.array([0.375, 0.375, 0.75 ]), 'K']]
    assert compare_path(stored_path,path)


