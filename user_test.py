from user_solver import Solver

# initiate variables
sigma_w = 1  # frequency bandwidth
omega_0 = 3  # central frequency
s = 10  # mesh points per wavelength
stability = 0.2  # time mesh stability factor

# initiate solver with user input variables
solver = Solver(s=s, stability=stability, simulation_time=1000)

# add pulse
solver.add_oscillating_pulse(sigma_w, (200, 150), omega_0, direction="right")

# add a reflecting square in top left
solver.add_reflect_square((3 * solver.size // 4, 3 * solver.size // 4), (-1, -1))

# change the boundary to have reflect on the bottom
solver.set_reflect_boundaries(up=True, down=True, left=True, right=True)

print(solver.steps)
print(solver.pulses[0].end_step())
# exit()

# add a material in the bottom right
# solver.add_material_square((3 * solver.size // 4, 3 * solver.size // 4), (4 * solver.size // 5, 4 * solver.size // 5),
#                            epsilon_rel=8)

solver.solve()
