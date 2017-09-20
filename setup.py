try:
      from setuptools import setup
except ImportError:
      from distutils.core import setup

setup(name='OpenDFT',
      version='1.0',
      py_modules=['main','solid_state_tools','exciting_handler'],
      author='Jannick Weisshaupt',
      author_email='jannickw@gmx.de',
    install_requires=['mayavi']
      )