from solver import Solver
from analyser import FileLoader
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

pulse = solver.add_oscillating_pulse(sigma_w, (1, 0.1), omega_0)
pulse.plot_time()

material = solver.create_material()

material.set_fixed_length_waveguide((1, 0.2), 5, 0.1, 0.5, 3)
material.plot()

solver.save("waveguide")
solver.solve(step_frequency=5)

loadfile = FileLoader("waveguide")
loadfile.play(1)
