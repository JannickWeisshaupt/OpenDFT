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

scipy_path = os.path.dirname(scipy.__file__) #use this if you are also using scipy in your application


def find_files(mypath,starts_with,exclude=()):

    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    vtk_files = []

    for file in onlyfiles:
        if file in exclude:
            continue
        if file.startswith(starts_with):
            vtk_files.append(join(mypath,file))
    return vtk_files

vtk_files = find_files('/usr/lib','libvtk')
hdf_files = find_files('/usr/lib/x86_64-linux-gnu','libhdf5_serial')
mpi1 = find_files('/usr/lib','libmpi')
mpi2 = find_files('/usr/lib/openmpi/lib','libmpi')

build_exe_options = {"packages": ["pyface.ui.qt4", "tvtk.vtk_module", "tvtk.pyface.ui.wx", "matplotlib.backends.backend_qt4",'pkg_resources._vendor','pkg_resources.extern','pygments.lexers',
                                  'tvtk.pyface.ui.qt4','pyface.qt','PySide','pyface.qt.QtGui','pyface.qt.QtCore','Shiboken'],
                     "include_files": [(str(scipy_path), "scipy"), #for scipy
                    (matplotlib.get_data_path(), "mpl-data"),'/home/jannick/python_programs/OpenDFT/data',
                                       '/usr/lib/x86_64-linux-gnu/libpyside-python2.7.so.1.2','/usr/lib/x86_64-linux-gnu/libpyside-python2.7.so.1.2.2',
                                       '/usr/lib/x86_64-linux-gnu/libshiboken-python2.7.so.1.2','/usr/lib/x86_64-linux-gnu/libshiboken-python2.7.so.1.2.2'
                                ,'/usr/lib/x86_64-linux-gnu/libpq.so','/usr/lib/x86_64-linux-gnu/libpq.so.5','/usr/lib/x86_64-linux-gnu/libpq.so.5.8',
                                       '/usr/lib/x86_64-linux-gnu/libmysqlclient.so','/usr/lib/x86_64-linux-gnu/libmysqlclient.so.20','/usr/lib/x86_64-linux-gnu/libmysqlclient.so.20.3.8',
                                       '/usr/lib/x86_64-linux-gnu/libnetcdf.so','/usr/lib/x86_64-linux-gnu/libnetcdf.so.11','/usr/lib/x86_64-linux-gnu/libnetcdf.so.11.0.0',
                                       '/usr/lib/x86_64-linux-gnu/libnetcdf_c++.so','/usr/lib/x86_64-linux-gnu/libnetcdf_c++.so.4','/usr/lib/x86_64-linux-gnu/libnetcdf_c++.so.4.2.0',
                                       '/usr/lib/libLSDyna.so.5.10','/usr/lib/libLSDyna.so.5.10.1','/usr/lib/x86_64-linux-gnu/libhdf5_serial_hl.so.10',
                                       '/usr/lib/x86_64-linux-gnu/libhdf5_serial_hl.so.10.0.2','/usr/lib/x86_64-linux-gnu/libsz.so','/usr/lib/x86_64-linux-gnu/libsz.so.2',
                                       '/usr/lib/x86_64-linux-gnu/libsz.so.2.0.1','/usr/lib/x86_64-linux-gnu/libaec.so','/usr/lib/x86_64-linux-gnu/libaec.so.0','/usr/lib/x86_64-linux-gnu/libaec.so.0.0.3',
                                       '/usr/lib/libgl2ps.so','/usr/lib/libgl2ps.so.0','/usr/lib/libgl2ps.so.0.0.0','/usr/lib/libVPIC.so.5.10','/usr/lib/libVPIC.so.5.10.1',
                                       '/usr/lib/libCosmo.so.5.10','/usr/lib/libCosmo.so.5.10.1','/usr/lib/libibverbs.so','/usr/lib/libibverbs.so.1','/usr/lib/libibverbs.so.1.0.0',
                                       '/usr/lib/libopen-rte.so.12','/usr/lib/openmpi/lib/libopen-pal.so.13.0.2','/usr/lib/openmpi/lib/libopen-pal.so','/usr/lib/libopen-pal.so.13',
                                       '/usr/lib/x86_64-linux-gnu/libhwloc.so.5','/usr/lib/x86_64-linux-gnu/libhwloc.so.5.6.8'
                                       ],

                     "includes":['numpy.core._methods', 'numpy.lib.format','PyQt4.QtCore','PyQt4.QtGui','pymatgen','pymatgen.symmetry.bandstructure'],
                     'excludes':'Tkinter',
                    "namespace_packages": ["ruamel.yaml"]
                    }

build_exe_options['include_files'].extend(vtk_files)
build_exe_options['include_files'].extend(hdf_files)
build_exe_options['include_files'].extend(mpi1)
build_exe_options['include_files'].extend(mpi2)


executables = [
    Executable('main.py', targetName="OpenDFT",icon="icon.ico")
]

setup(name='OpenDFT',
      version='1.0',
      description='OpenDFT',
      options = {"build_exe": build_exe_options},
      executables=executables,
      )