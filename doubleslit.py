from solver import Solver
from analyser import FileLoader
import matplotlib.pyplot as plt
import math

# initiate variables
sigma_w = 5.0 * 10 ** 9  # frequency bandwidth
omega_0 = 2 * 10 ** 9  # central frequency
s = 10  # mesh points per wavelength
stability = 1 / 10  # time mesh stability factor
C = 2.99 * (10 ** 8)

wavelength = 2 * math.pi * C / omega_0  # c = f * lambda, f = omega / 2pi

simulation_size = 4
simulation_time = 1.5 * simulation_size / C

print(simulation_time)

# initiate solver with user input variables
solver = Solver(points_per_wavelength=s, stability=stability, eps_r_max=1, mu_r_max=1, simulation_size=simulation_size, simulation_time= simulation_time)

pulse1 = solver.add_oscillating_pulse(sigma_w, (2, 0.1), omega_0, direction="right")

# pulse1.plot()
pulse1.plot_frequency()
# pulse1.plot_frequency_fft()
exit()

material = solver.create_material()

# double slits where slit separation is wavelength, slit width is wavelength
slit_space_upper = (simulation_size / 2) - wavelength / 2
slit_space_lower = (simulation_size / 2) + wavelength / 2

solver.set_reflect_square((0, 1), (1.62, 1.05))
solver.set_reflect_square((1.76, 1), (2.24, 1.05))
solver.set_reflect_square((2.38, 1), (simulation_size, 1.05))
solver.set_reflect_boundaries(up=False, down=False, left=False, right=False)

solver.save("far_spaced_double_slit")
solver.solve(realtime=False)

fileLoad = FileLoader("far_spaced_double_slit")
fileLoad.play(5)
