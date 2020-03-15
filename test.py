# Used for testing new features

from solver import Solver
from analyser import FileLoader
import math
import json

# initiate variables
sigma_w = 1 * 10 ** 9  # frequency bandwidth
omega_0 = 5 * 10 ** 9  # central frequency

print(2 * math.pi * 2.99 * (10 ** 8) / omega_0)
s = 10  # mesh points per wavelength
stability = 0.2  # time mesh stability factor

# initiate solver with user input variables
solver = Solver(points_per_wavelength=s, stability=stability, eps_r_max=4, mu_r_max=1, simulation_size=3,
                simulation_time=2.5 * 10 ** (-8))

# add pulses (a pulse must be added to run the simulation)
solver.add_oscillating_pulse(sigma_w, (0.8, 0.4), omega_0)  # this adds a point pulse

mat = solver.create_material()

mat.set_material_rect((2, 2), (2.8, 2.8), 3)

solver.save('test_json')
solver.solve(realtime=False)

fileLoad = FileLoader('test_json')
fileLoad.play()

with open ('test_json.txt') as json_file:
    data = json.load(json_file)

print("Number of simulation steps: {}".format(int(data['end_time'] / (data['dt'] * data['step_frequency']))))
