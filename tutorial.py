from solver import Solver, load_file

# initiate variables
sigma_w = 1 * 10 ** 9  # frequency bandwidth
omega_0 = 5 * 10 ** 9  # central frequency

s = 10  # mesh points per wavelength
stability = 0.2  # time mesh stability factor

# initiate solver with user input variables
solver = Solver(points_per_wavelength=s, stability=stability, eps_r_max=4, mu_r_max=1, simulation_size=3, simulation_time=2.5 * 10**(-8))

# add pulses (a pulse must be added to run the simulation)
solver.add_oscillating_pulse(sigma_w, (0.8, 0.4), omega_0)  # this adds a point pulse

# REQUIRED
material = solver.create_material()

# OPTIONAL
solver.set_reflect_square((0.98, 0.8), (solver.size, 1))

# OPTIONAL
solver.add_data_collector(0.82, 0.95,)
solver.add_data_collector(1.9, 2.13)

# OPTIONAL
material.set_fixed_length_waveguide(wire_length=5, thickness=0.1, start_point=(0.2, 1.5), curved_ratio=1, epsilon_rel=4)
material.plot()

# OPTIONAL
# change the boundary to be reflective
solver.set_reflect_boundaries(up=False, down=False, left=False, right=False)

# If you wish to save as an npy file, uncomment below and use solver.load() at the end to see simulation (recommended)
# solver.save('test')

# If you wish to see the simulation...
solver.solve()

# Use this if you have used solver.save
# solver.load()

# Work in progress
# for collector in solver.data_collectors:
#     collector.plot()
#     collector.fft()
