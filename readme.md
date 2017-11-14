# OpenDFT

OpenDFT is a free and open source software that brings cutting edge solid state research to the people. It is a graphical program that interacts with various scientific terminal based solid state software packages. It visualizes inputs, such as the crystal structure and visualizes outputs, such as band structures and optical spectra. OpenDFT is designed to be DFT engine agnostic so that you can easily switch between different scientific codes and compare results seamlessly.

### Installing

OpenDFT runs on windows and linux and python versions 2 and 3. In order to be able to run it you need to install its dependencies.

Ubuntu, Linux Mint, etc. (linux distros with aptitude package manager):

sudo apt-get install python-numpy <br>
sudo apt-get install python-pip <br>
sudo apt-get install python-scipy <br>
sudo apt-get install python-matplotlib <br>
sudo apt-get install python-dev <br>
sudo apt-get install python-qt4 <br>
sudo apt-get install python-vtk python-wxgtk2.6 python-setuptools python-configobj <br>
pip install mayavi (falls das nicht geht: "sudo pip install mayavi" und im folgenden immer sudo pip) <br>
pip install pyface <br>
pip install six <br>
pip install traitsui <br>
pip install pygments <br>
pip install periodictable <br>
sudo apt-get install git <br>

sudo apt-get install nwchem <br>
sudo apt-get install quantum-espresso <br>

cd ordner/in/den/es/installiert/werden/soll <br>
git clone https://github.com/JannickWeisshaupt/OpenDFT.git <br>
cd OpenDFT <br>
./main.py <br>


## Getting Started