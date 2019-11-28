from user_solver import Solver

# initiate variables
sigma_w = 0.5  # frequency bandwidth
omega_0 = 5  # central frequency
s = 10  # mesh points per wavelength
stability = 0.1  # time mesh stability factor

# initiate solver with user input variables
solver = Solver(s=s, stability=stability, eps_r_max=5, mu_r_max=1, simulation_time=1000)
# solver.add_oscillating_pulse(sigma_w, (100, 50), omega_0, direction="right")
solver.add_oscillating_pulse(sigma_w, (100, 50), omega_0)

# add a reflecting square in top left
solver.add_reflect_square((3 * solver.size // 4, 3 * solver.size // 4), (-1, -1))

solver.add_material_convex((100, 100), 100, 80, epsilon_rel=5, mu_rel=1)

# add a material in the bottom right
# solver.add_material_square((3 * solver.size // 4, 3 * solver.size // 4), (4 * solver.size // 5, 4 * solver.size // 5),
#                            epsilon_rel=8)

# change the boundary to have reflect
solver.set_reflect_boundaries(up=True, down=True, left=False, right=False)

solver.solve()
