from user_solver import Solver

# initiate variables
sigma_w = 1  # frequency bandwidth
omega_0 = 3  # central frequency
s = 10  # mesh points per wavelength
stability = 0.1  # time mesh stability factor

# initiate solver with user input variables
solver = Solver(points_per_wavelength=s, stability=stability, eps_r_max=4, mu_r_max=1, simulation_time=1000)

# add pulses (a pulse must be added to run the simulation)
solver.add_oscillating_pulse(sigma_w, (100, 50), omega_0, direction="right")  # this adds a point pulse

# add a reflecting square in bottom
solver.add_reflect_square((3 * solver.size // 4, 2 * solver.size // 4), (-1, 3 * solver.size // 4))

# solver.set_material_convex((100, 100), 100, 80, epsilon_rel=5, mu_rel=1)
solver.set_waveguide((solver.size//2, 3 * solver.size//4), 50, 50, 30, epsilon_rel=4)
solver.add_reflect_square((0, 174), (100, 175))
solver.add_reflect_square((130, 174), (-1, 175))

# add a material in the top right
solver.set_material_rect((0, 3 * solver.size // 4), (1 * solver.size // 5, solver.size),
                         epsilon_rel=4)

solver.plot_materials()

solver.plot_materials()
# change the boundary to have reflect
# solver.set_reflect_boundaries(up=False, down=False, left=False, right=False)

# run
solver.solve()
