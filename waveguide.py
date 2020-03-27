from solver import Solver
import math

# initiate variables
sigma_w = 1.0 * 10 ** 9  # frequency bandwidth
omega_0 = 6 * 10 ** 9  # central frequency
s = 10  # mesh points per wavelength
stability = 1 / 10  # time mesh stability factor
C = 2.99 * (10 ** 8)

wavelength = 2 * math.pi * C / omega_0  # c = f * lambda, f = omega / 2pi

simulation_size = 4
simulation_time = 5 * simulation_size / C

print("Simulation time (seconds):", simulation_time)

# initiate solver with user input variables
solver = Solver(points_per_wavelength=s, stability=stability, eps_r_max=3, mu_r_max=1, simulation_size=simulation_size, simulation_time= simulation_time)

solver.set_reflect_boundaries(False, False, False, False)

pulse = solver.add_oscillating_pulse(sigma_w, (0.7, 0.2), omega_0)

material = solver.create_material()
material.set_fixed_length_waveguide((0.5, 0.25), 4, 0.7, 0.4, 3)
material.set_fixed_length_waveguide((0.65, 0.25), 4, 0.7, 0.1, 1)

solver.set_reflect_square((0, 0.25), (0.5, 0.35))
solver.set_reflect_square((0.9, 0.25), (4, 0.35))

solver.save("waveguide")
solver.solve(step_frequency=5)
