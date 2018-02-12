import sys
import os
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

scipy_path = os.path.dirname(scipy.__file__) #use this if you are also using scipy in your application

build_exe_options = {"packages": ["pyface.ui.qt4", "tvtk.vtk_module", "tvtk.pyface.ui.wx", "matplotlib.backends.backend_qt4",'pkg_resources._vendor','pkg_resources.extern','pygments.lexers',
                                  'tvtk.pyface.ui.qt4','pyface.qt','PySide','pyface.qt.QtGui','pyface.qt.QtCore','Shiboken'],
                     "include_files": [(str(scipy_path), "scipy"), #for scipy
                    (matplotlib.get_data_path(), "mpl-data"),'/home/jannick/python_programs/OpenDFT/data',
                                       '/usr/lib/x86_64-linux-gnu/libpyside-python2.7.so.1.2','/usr/lib/x86_64-linux-gnu/libpyside-python2.7.so.1.2.2',
                                       '/usr/lib/x86_64-linux-gnu/libshiboken-python2.7.so.1.2','/usr/lib/x86_64-linux-gnu/libshiboken-python2.7.so.1.2.2',
                                       '/usr/lib/libvtkCommonPythonD.so.5.10','/usr/lib/libvtkCommonPythonD.so.5.10.1'],
                     "includes":['numpy.core._methods', 'numpy.lib.format','PyQt4.QtCore','PyQt4.QtGui','pymatgen','pymatgen.symmetry.bandstructure'],
                     'excludes':'Tkinter',
                    "namespace_packages": ["ruamel.yaml"]
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