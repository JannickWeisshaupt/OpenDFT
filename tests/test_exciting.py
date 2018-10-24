import sys, os
import pytest

myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, myPath + '/../')

import numpy as np
import src.solid_state_tools as sst
from src.exciting_handler import Handler

unit_cell = 6.719 * np.array([[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]])
atoms = np.array([[0.0, 0.0, 0.0, 6], [0.25, 0.25, 0.25, 6]])
diamond = sst.CrystalStructure(unit_cell, atoms)

if not os.path.isdir(os.path.dirname(os.path.abspath(__file__))+'/outputs/'):
    os.mkdir(os.path.dirname(os.path.abspath(__file__))+'/outputs/')

def start_engine(blocking=True):
    pass

handler = Handler()
handler._start_engine = start_engine
handler.project_directory = os.path.dirname(os.path.abspath(__file__))+'/outputs/'

def test_scf():
    handler.start_ground_state(diamond)
    reread_diamond = handler.parse_input_file(handler.project_directory+'/exciting_files/input.xml')

    all_tests = []

    dens = reread_diamond.density()

    dens_test = np.abs(dens - 0.31677024692785494) < 1e-12
    all_tests.append(dens_test)

    info = reread_diamond.lattice_information()
    diamond_info = {'crystal system': 'cubic', 'point group': 'm-3m', 'space group': 'Fd-3m'}

    info_test = info == diamond_info
    all_tests.append(info_test)

    assert all(all_tests)

