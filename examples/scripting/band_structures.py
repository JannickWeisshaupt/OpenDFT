
engine.start_ground_state(structure, blocking=True)  # start the calculation
band_structure = engine.read_bandstructure()
add_data(band_structure,'test')