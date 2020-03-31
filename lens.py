from solver import Solver
from analyser import FileLoader

# initiate variables
sigma_w = 1 * 10 ** 9  # frequency bandwidth
omega_0 = 5 * 10 ** 9  # central frequency
s = 10  # mesh points per wavelength
stability = 0.1  # time mesh stability factor

# initiate solver with user input variables
solver = Solver(points_per_wavelength=s, stability=stability, eps_r_max=5, mu_r_max=1, simulation_size=3, simulation_time=2 * 10**(-8))
pulse = solver.add_sinusoidal_pulse(location=(1.5, 0.2), omega_0=omega_0, direction="right")
material = solver.create_material()
material.set_material_convex((1.5, 1), radius=2.2, thickness=1.2, epsilon_rel=5)

material.plot()

solver.set_reflect_boundaries(up=False, down=False, right=False, left=False)

solver.save('lens')
solver.solve()

fileloader = FileLoader('lens')
fileloader.play(1)
