try:
      from setuptools import setup
except ImportError:
      from distutils.core import setup

setup(name='opendft',
      version='1.0',
      py_modules=['main','solid_state_tools','exciting_handler','abinit_handler','quantum_espresso_handler','nwchem_handler','syntax','TerminalClass','visualization','little_helpers'],
      author='Jannick Weisshaupt',
      author_email='jannickw@gmx.de',
      install_requires=['setuptools','numpy','matplotlib','periodictable','pyface','six','pymatgen','PySide','mayavi'],
      maintainer='Jannick Weisshaupt',
      maintainer_email='jannickw@gmx.de',
      description='A Gui application for density functional theory calculations'
      )