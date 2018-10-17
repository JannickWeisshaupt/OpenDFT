# Welcome to the OpenDFT scripting console
#
# You can normally(*) use python to script here in addition to OpenDFT objects that can be
# used to calculate electronic properties and visualize them. The structure, engine options etc. are loaded, when the
# scripting window is opened and can be changed within the script.
#
# Predefined variables:
#
# engine:               Dft engine object that can be used to calculate electronic properties.
#
# structure:            crystal or molecular structure from the main application.
#                       If you defined a structure with the main window you can directly use it here.
#
# plot_structure:       Function that updates the plot in the main window.
#
# plot_scf:             Function that plots the current scf convergence in the main window.
#
# CrystalStructure:     Class to describe a crystal structure.
#
# MolecularStructure:   Class to describe a molecule
#
# Use the help function to learn more about the variables,
# e.g. help(engine) and help(engine.start_ground_state) should be quite helpful
#
# (*) For technical reasons matplotlib can be used but has to be put at the end of the script after a special seperator,
# namely [dollarsign]matplotlib (see example below). This seperator must be used only once. All code after it will run in the
# main thread, i.e. the application freezes while it runs. Make sure the computationally heavy stuff comes before the seperator.
#
# The following code is an runnable (press f5) example of how to convergence a exciting calculation with
# respect to the k vector sampling and rgkmax (with very few steps for faster run-time)
#
import numpy as np


k_vecs = range(4,8)
basis_parameters = np.linspace(4,8,6)

energies = np.zeros((len(k_vecs),len(basis_parameters)))

engine.scf_options['do'] = 'fromfile'

for i,k_vec in enumerate(k_vecs):
    for j,basis in enumerate(basis_parameters):
        # set the options
        engine.scf_options['rgkmax'] = '{0:1.3f}'.format(basis)
        engine.scf_options['ngridk'] = '{0} {0} {0}'.format(k_vec)

        engine.start_ground_state(structure, blocking=True)  # start the calculation

        try:
            scf_list = engine.read_scf_status()  # Read the full list of scf energies
            plot_scf(scf_list)
            energy = scf_list[-1, 1]  # append only the last to the energies list
            energies[i,j] = energy
        except Exception as e:
            print(repr(e))
            energies[i, j] = np.nan

np.savetxt(engine.project_directory+'/convergence.dat',energies)

## Plot with matplotlib

$matplotlib

xv, yv = np.meshgrid(basis_parameters,k_vecs)


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
plt.ion()

fig = plt.figure(1,figsize=(10,8))
fig.clf()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(xv,yv,energies-energies.min(),cmap='viridis_r')
ax.set_ylabel('k sampling')
ax.set_xlabel('rgkmax')
plt.tight_layout()
plt.savefig(engine.project_directory+'/convergence_3d.pdf')

fig2 = plt.figure(2,figsize=(10,8))
fig2.clf()
ax1 = fig2.add_subplot(211)
ax2 = fig2.add_subplot(212)
ax1.semilogy(k_vecs,energies[:,-1]-energies.min(),'-o')
ax1.set_xlabel('k sampling')
ax2.semilogy(basis_parameters,energies[-1,:]-energies.min(),'-o')
ax2.set_xlabel('rgkmax')
plt.tight_layout()
plt.savefig(engine.project_directory+'/convergence.pdf')

plt.show()
