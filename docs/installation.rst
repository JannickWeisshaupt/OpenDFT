#######################
How to install OpenDFT
#######################
OpenDFT is very easy to install on its own. Just follow the steps below and it should work. OpenDFT relies on DFT codes, e.g. exciting, abinit etc.
that are already installed on your system. They must be available in the terminal that start OpenDFT with their standard name, i.e. **excitingser** for exciting, abinit for abinit etc.,
otherwise they cannot be used.

- Download the latest release from the OpenDFT release_ page and save it on your work station.
- Extract the archive (for linux use the .tar.gz one) somewhere on your machine
- Add the OpenDFT executable in the OpenDFT directory to your PATH variable or make an alias, e.g for bash open the ~/.bashrc file and add:|br| ``alias OpenDFT=/home/user/path/to/OpenDFT/OpenDFT``
- Reload the environmental variable of the terminal with ``source ~/.bashrc`` or simply restart the terminal.
- Type ``OpenDFT`` in your terminal to start OpenDFT

If OpenDFT does not start and some errors appear in the terminal try to install qt4 on your system, for e.g. ubuntu use |br| **sudo apt install qt4-default**.


#######################
DFT codes installation
#######################
OpenDFT is a graphical user interface to command line quantum chemistry and/or DFT codes. As such it requires those to be installed and accessible through their standard
names in the terminal that started OpenDFT, i.e. through the environmental variables supplied to OpenDFT. Below is a "how to install" list of the supported codes.

=======================
Nwchem
=======================
``sudo apt-get install nwchem``
or download the nwchem-release_ and follow these instructions_

=======================
Quantum espresso
=======================
``sudo apt-get install quantum-espresso``

.. _release: https://github.com/JannickWeisshaupt/OpenDFT/releases
.. _nwchem-release: https://github.com/nwchemgit/nwchem/releases/download/v6.8-release/nwchem-6.8-release.revision-v6.8-47-gdf6c956-src.2017-12-14.tar.bz2
.. _instructions: http://www.nwchem-sw.org/index.php/Compiling_NWChem
.. |br| raw:: html

    <br />