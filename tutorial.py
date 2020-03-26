from solver import Solver
from analyser import FileLoader

# initiate variables
sigma_w = 1 * 10 ** 9  # frequency bandwidth
omega_0 = 5 * 10 ** 9  # central frequency

s = 10  # mesh points per wavelength
stability = 0.2  # time mesh stability factor

# initiate solver with user input variables
solver = Solver(points_per_wavelength=s, stability=stability, eps_r_max=4, mu_r_max=1, simulation_size=3, simulation_time=2.5 * 10**(-8))

# REQUIRED - add a pulse
my_pulse = solver.add_oscillating_pulse(sigma_w, (0.8, 0.4), omega_0)  # this adds a point pulse
my_pulse.plot_time()

# REQUIRED - create the background material
my_material = solver.create_material()

# OPTIONAL - creating reflecting shapes
solver.set_reflect_square((0.98, 0.8), (solver.size, 1))

# OPTIONAL - create the material
my_material.set_fixed_length_waveguide(wire_length=5, thickness=0.1, start_point=(0.2, 1.5), curved_ratio=1, epsilon_rel=4)
my_material.plot_time()

# OPTIONAL - change the boundary to be reflective
solver.set_reflect_boundaries(up=False, down=False, left=False, right=False)

# OPTIONAL - save file (REQUIRED for data analysis)
save_file_name = 'test'
solver.save(save_file_name)

# If you wish to see the simulation...
solver.solve()

# Use this if you have used solver.save
file_loader = FileLoader(save_file_name)
file_loader.play(interval=50)

# Creates data collector object in the middle of simulation
data1 = file_loader.create_data_collector((file_loader.length_y / 2, file_loader.length_x / 2))
data1.plot_frequency()
data1.plot_time()
