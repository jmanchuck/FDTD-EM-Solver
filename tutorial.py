from user_solver import Solver

# initiate variables
sigma_w = 1 * 10**(9)  # frequency bandwidth
omega_0 = 5 * 10**(9)  # central frequency
s = 10  # mesh points per wavelength
stability = 0.1  # time mesh stability factor

# initiate solver with user input variables
solver = Solver(points_per_wavelength=s, stability=stability, eps_r_max=4, mu_r_max=1, simulation_size=3, simulation_time=1 * 10**(-8))

# add pulses (a pulse must be added to run the simulation)
# solver.add_oscillating_pulse(sigma_w, (100, 50), omega_0, direction="right")  # this adds a plane pulse
solver.add_oscillating_pulse(sigma_w, (0.6, 0.6), omega_0)  # this adds a point pulse

# solver.set_material_convex((100, 100), 100, 80, epsilon_rel=5, mu_rel=1)
# solver.set_waveguide((solver.size//2, solver.size//2), radius=80, rectangle_length=130, thickness=30, epsilon_rel=4)

# add a material in the top left
solver.set_material_rect((0, 0), (0.5, 0.5), epsilon_rel=4)

solver.plot_materials()
# change the boundary to have reflect
# solver.set_reflect_boundaries(up=False, down=False, left=False, right=False)

# run
solver.solve(save=True)
