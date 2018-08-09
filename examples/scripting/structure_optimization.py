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
# The following code is an runnable (press f5) example of how to find the optimal cell scale (with very few steps for faster run-time)
#
import numpy as np

atoms = structure.atoms  # Save the current atom information for future use
lattice_vectors = structure.lattice_vectors  # save the lattice vectors

energies = []
scales = np.linspace(0.8, 1.2, 5)  # scale factor for unit cell
scales_succes = []

for scale in scales:
    structure_temp = CrystalStructure(scale * lattice_vectors, atoms)
    plot_structure(structure_temp)  # plot the structure in the main window
    engine.start_ground_state(structure_temp, blocking=True)  # start the calculation
    try:
        scf_list = engine.read_scf_status()  # Read the full list of scf energies
        plot_scf(scf_list)
        energies.append(scf_list[-1, 1])  # append only the last to the energies list
        scales_succes.append(scale)
    except Exception as e:
        print(repr(e))

print('Emin = {0:1.2f}'.format(min(energies)))
optimal_scale = scales_succes[energies.index(min(energies))]
print('Optimal cell scale = {0:1.2f}'.format(optimal_scale))

plot_structure(structure)

## Plot with matplotlib

$matplotlib

import matplotlib.pyplot as plt
plt.ion()

plt.figure(1)
plt.clf()
plt.plot(scales_succes, energies, 'k', linewidth=2)
plt.show()