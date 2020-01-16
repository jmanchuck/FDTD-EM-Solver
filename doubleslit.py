from user_solver import Solver

# initiate variables
sigma_w = 1 * 10 ** 9  # frequency bandwidth
omega_0 = 20 * 10 ** 9  # central frequency
s = 10  # mesh points per wavelength
stability = 0.1  # time mesh stability factor

# initiate solver with user input variables
solver = Solver(points_per_wavelength=s, stability=stability, eps_r_max=4, mu_r_max=1, simulation_size=2, simulation_time=1.3 * 10**(-8))

solver.add_oscillating_pulse(sigma_w, (2.5, 0.1), omega_0, direction="right")

solver.set_reflect_square((0, 0.5), (0.85, 0.51))
solver.set_reflect_square((0.95, 0.5), (1.05, 0.51))
solver.set_reflect_square((1.15, 0.5), (5, 0.51))

solver.solve(realtime=True)
