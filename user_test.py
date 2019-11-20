from user_solver import Solver

# initiate variables
sigma_w = 1  # frequency bandwidth
omega_0 = 5  # central frequency
s = 10  # mesh points per wavelength
stability = 0.2  # time mesh stability factor

# initiate solver with user input variables
solver = Solver(s=s, stability=stability, simulation_time=1000)

# add pulse
solver.add_oscillating_pulse(sigma_w, (150, 150), omega_0)

# add a reflecting square in top left
solver.add_reflect_square((0, 0), (solver.size // 4, solver.size // 4))

# change the boundary to have reflect on the bottom
solver.set_reflect_boundaries(up=False, down=True, left=False, right=False)

# add a material in the bottom right
solver.add_material_square((3 * solver.size // 4, 3 * solver.size // 4), (4 * solver.size // 5, 4 * solver.size // 5),
                           epsilon_rel=8)

solver.solve()
