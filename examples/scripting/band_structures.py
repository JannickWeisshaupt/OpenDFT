import numpy as np

basis_parameters = np.linspace(4,8,3)

for basis_parameter in basis_parameters:
    engine.scf_options['rgkmax'] = '{0:1.3f}'.format(basis_parameter)
    engine.start_ground_state(structure, blocking=True,band_structure_points=k_path)  # start the calculation
    band_structure = engine.read_bandstructure()
    add_data(band_structure,'rgkmax_{0:1.2f}'.format(basis_parameter))