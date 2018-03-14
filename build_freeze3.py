import sys
import os
from os import listdir
from os.path import isfile, join
from cx_Freeze import setup, Executable
import cx_Freeze.hooks
def hack(finder, module):
    return
cx_Freeze.hooks.load_matplotlib = hack
import scipy
import matplotlib
from encodings import ascii
from encodings import idna
from encodings import unicode_escape
import ruamel.yaml



scipy_path = os.path.dirname(scipy.__file__) #use this if you are also using scipy in your application





build_exe_options = {"packages": ["pyface.ui.qt4", "tvtk.vtk_module", "tvtk.pyface.ui.wx", "matplotlib.backends.backend_qt4",'pkg_resources._vendor','pkg_resources.extern','pygments.lexers',
                                  'tvtk.pyface.ui.qt4','pyface.qt','pyface.qt.QtGui','pyface.qt.QtCore','numpy','matplotlib','mayavi'],
                     "include_files": [(str(scipy_path), "scipy"), #for scipy
                    (matplotlib.get_data_path(), "mpl-data"),'/home/jannick/python_programs/OpenDFT/data','/home/jannick/.local/lib/python3.4/site-packages/pymatgen/core/bond_lengths.json',
                                       '/home/jannick/.local/lib/python3.4/site-packages/pymatgen/core/func_groups.json','/home/jannick/.local/lib/python3.4/site-packages/pymatgen/core/libxc_docs.json',
                                       '/home/jannick/.local/lib/python3.4/site-packages/pymatgen/core/periodic_table.json','/home/jannick/.local/lib/python3.4/site-packages/pymatgen/core/reconstructions_archive.json',
                                       '/usr/local/lib/python3.4/site-packages/PyQt4','/home/jannick/.local/lib/python3.4/site-packages/mayavi','/home/jannick/.local/lib/python3.4/site-packages/ruamel',
                                       '/home/jannick/.local/lib/python3.4/site-packages/pyface','/home/jannick/.local/lib/python3.4/site-packages/tvtk',
                                       '/home/jannick/.local/lib/python3.4/site-packages/pymatgen/core/bond_lengths.json','/home/jannick/.local/lib/python3.4/site-packages/pymatgen/core/func_groups.json',
                                       '/home/jannick/.local/lib/python3.4/site-packages/pymatgen/core/libxc_docs.json','/home/jannick/.local/lib/python3.4/site-packages/pymatgen/core/periodic_table.json',
                                       '/home/jannick/.local/lib/python3.4/site-packages/pymatgen/core/reconstructions_archive.json'
                                       ],

                     "includes":['PyQt4.QtCore','PyQt4.QtGui','pymatgen','pymatgen.symmetry.bandstructure','mayavi','PyQt4'],
                     'excludes':'Tkinter',
                    "namespace_packages": ['mayavi']
                    }


executables = [
    Executable('main.py', targetName="OpenDFT",icon="icon.ico")
]

setup(name='OpenDFT',
      version='1.0',
      description='OpenDFT',
      options = {"build_exe": build_exe_options},
      executables=executables,
      )