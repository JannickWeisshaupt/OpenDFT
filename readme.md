OpenDFT is a free and open source software that brings cutting edge solid state research to the people. It is a graphical program that interacts with various scientific terminal based solid state software packages. It visualizes inputs, such as the crystal structure and visualizes outputs, such as band structures and optical spectra. OpenDFT is designed to be DFT engine agnostic so that you can easily switch between different scientific codes and compare results seamlessly.

Installation:

OpenDFT runs on windows and linux and python versions 2 and 3. In order to be able to run it you need to install its dependencies.

Ubuntu, Linux Mint, etc. (linux distros with aptitude package manager):

sudo apt-get install python-numpy
sudo apt-get install python-pip
sudo apt-get install python-scipy
sudo apt-get install python-matplotlib
sudo apt-get install python-dev
sudo apt-get install python-qt4
sudo apt-get install python-vtk python-wxgtk2.6 python-setuptools python-configobj
pip install mayavi (falls das nicht geht: "sudo pip install mayavi" und im folgenden immer sudo pip)
pip install pyface
pip install six
pip install traitsui
pip install pygments
pip install periodictable
sudo apt-get install git

sudo apt-get install nwchem
sudo apt-get install quantum-espresso

cd ordner/in/den/es/installiert/werden/soll
git clone https://github.com/JannickWeisshaupt/OpenDFT.git
cd OpenDFT
./main.py 