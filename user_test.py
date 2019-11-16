from user_solver import Solver

# initiate variables
sigma_w = 1  # frequency bandwidth
omega_0 = 0  # central frequency
s = 10  # mesh points per wavelength
stability = 0.2  # time mesh stability factor

# initiate solver with user input variables
solver = Solver(sigma_w=sigma_w, omega_0=omega_0, s=s, stability=stability)
print(solver.size)

# add a reflecting square in top left
solver.add_reflect_square((0, 0), (solver.size // 4, solver.size // 4))

# change the boundary to have reflect on the bottom
solver.set_boundaries(up=False, down=True, left=False, right=False)

solver.solve()
