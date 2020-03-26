from solver import Solver
from analyser import FileLoader
import math

# initiate variables
sigma_w = 3.0 * 10 ** 9  # frequency bandwidth
omega_0 = 6 * 10 ** 9  # central frequency
s = 10  # mesh points per wavelength
stability = 1 / 10  # time mesh stability factor
C = 2.99 * (10 ** 8)

wavelength = 2 * math.pi * C / omega_0  # c = f * lambda, f = omega / 2pi

simulation_size = 4
simulation_time = 1.5 * simulation_size / C

print("Simulation time (seconds):", simulation_time)

# initiate solver with user input variables
solver = Solver(points_per_wavelength=s, stability=stability, eps_r_max=1, mu_r_max=1, simulation_size=simulation_size, simulation_time= simulation_time)

pulse1 = solver.add_oscillating_pulse(sigma_w, (2, 0.1), omega_0, direction="right")

material = solver.create_material()

# double slits where slit separation is wavelength, slit width is wavelength
slit_space_upper = (simulation_size / 2) - wavelength / 2
slit_space_lower = (simulation_size / 2) + wavelength / 2

solver.set_reflect_square((0, 1), (1.62, 1.05))
solver.set_reflect_square((1.76, 1), (2.24, 1.05))
solver.set_reflect_square((2.38, 1), (simulation_size, 1.05))
solver.set_reflect_boundaries(up=False, down=False, left=False, right=False)

solver.save("d_slit_test")
solver.solve(realtime=False, step_frequency=2)

fileLoad = FileLoader("d_slit_test")
fileLoad.play(1)
