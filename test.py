# Used for testing new features

from solver import Solver
from analyser import FileLoader

sigma_w = 10 * 10 ** 9  # frequency bandwidth
omega_0 = 10 * 10 ** 9  # central frequency

print(2.99 * (10 ** 8) / omega_0)
s = 10  # mesh points per wavelength
stability = 0.2  # time mesh stability factor

solver = Solver(points_per_wavelength=s, stability=stability, eps_r_max=4, mu_r_max=1, simulation_size=3,
                simulation_time=2.5 * 10 ** (-8))

solver.add_oscillating_pulse(sigma_w, (0.1, 0.05), omega_0)

for pulse in solver.pulses:
    pulse.plot()
# mat = solver.create_material()
# mat.plot()
#
# solver.save("test")
# solver.solve(realtime=False)
#
# fileloader = FileLoader('test')
# fileloader.play()