from solver import Solver

# initiate variables
sigma_w = 1 * 10 ** 9  # frequency bandwidth
omega_0 = 5 * 10 ** 9  # central frequency
s = 10  # mesh points per wavelength
stability = 0.3  # time mesh stability factor

# initiate solver with user input variables
solver = Solver(points_per_wavelength=s, stability=stability, eps_r_max=5, mu_r_max=1, simulation_size=3, simulation_time=2.5 * 10**(-8))
solver.add_oscillating_pulse(sigma_w, (1.5, 0.1), omega_0, direction="right")
material = solver.create_material()
material.set_material_convex((1.5, 1.5), radius=2, thickness=1, epsilon_rel=5)
solver.set_reflect_boundaries(up=True, down=True)

solver.save('lens')
solver.solve(realtime=False)
solver.load()
