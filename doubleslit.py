from user_solver import Solver

# initiate variables
sigma_w = 2  # frequency bandwidth
omega_0 = 30  # central frequency
s = 10  # mesh points per wavelength
stability = 0.1  # time mesh stability factor

# initiate solver with user input variables
solver = Solver(points_per_wavelength=s, stability=stability, eps_r_max=1, mu_r_max=1, simulation_size=1400000000, simulation_time=1000)

solver.add_oscillating_pulse(sigma_w, (150, 20), omega_0, direction="right")

solver.add_reflect_square((0, 40), (120, 42))
solver.add_reflect_square((140, 40), (160, 42))
solver.add_reflect_square((180, 40), (solver.size, 42))

solver.solve()
