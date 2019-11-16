import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd
from time import process_time

# Physical constants
c = 2.99 * (10 ** 8)
mu = 1.26 * (10 ** (-6))
epsilon = 8.85 * (10 ** (-12))

# ALL USER INPUTS HERE
sigma_w = 1         # frequency bandwidth
omega_0 = 0         # central frequency
s = 10              # mesh points per wavelength
stability = 0.2     # time mesh stability factor

epsilon_relative = 2  # change this to a list in the future
mu_relative = 1       # change this to a list

# improvements:
# - take into account epsilon max to determine how we will sample a wavelength (to get dx or ds)
# - user inputs real time simulation length in seconds, use dt to calculate how many steps


# Derived constants, calculate from user input but stays constant throughout simulation
n_max = math.sqrt(epsilon_relative * mu_relative)
sigma_t = 1 / sigma_w
omega_max = omega_0 + 3 * sigma_w
lambda_min = math.pi * 2 * c / (n_max * omega_max)
l_x, l_y = 30 * lambda_min, 30 * lambda_min     # size of simulation is 50 wavelengths
ds = lambda_min / s
dt = min(ds * stability / c, (4 * math.pi) / omega_max)

# Simulation time
steps = 1500
end_time = dt * steps

# Gaussian derivative properties
pulseMid = 3 * sigma_t  # the middle of the gaussian derivative (approximate)


def plot_gaussian_derivative():

    data = []
    t = 0

    while t < 5 * pulseMid:
        gauss = (-t + pulseMid) * (1 / (sigma_t * math.sqrt(2 * math.pi))) * (
            math.exp(-((t - pulseMid) ** 2) / (2 * (sigma_t ** 2))))
        t += dt
        data.append(gauss)
    print(len(data))
    plt.plot(data)
    plt.show()
    exit()


def gaussian_derivative_maxima():
    t = 0
    current = 0
    prev = 0
    while t < pulseMid:
        current = (-t + pulseMid) * (1 / (sigma_t * math.sqrt(2 * math.pi))) * (
            math.exp(-((t - pulseMid) ** 2) / (2 * (sigma_t ** 2))))
        if current < prev:
            return prev
        t += dt
        prev = current
    return current


maxima = gaussian_derivative_maxima()

# Simulation variables, these change throughout the simulation
time = 0
step = 0
size = int(l_x / ds)
h = np.zeros((size, size))
ex = np.zeros((size + 1, size))
ey = np.zeros((size, size + 1))

# animation plotting properties
frames = np.arange(0, end_time, dt)
fig, ax = plt.subplots()

# material properties matrix
eps_arr = (dt / ds) * (np.ones((size, size)) / epsilon)
mu_arr = (dt / ds) * (np.ones((size, size)) / mu)

# adding square material in top right
eps_arr[:(size // 3), 2 * (size // 3):] *= epsilon_relative
mu_arr[:(size // 3), 2 * (size // 3):] *= 1

# reflecting box in the bottom left
row_upper = 2 * (size // 3)
row_lower = size - 1
col_left = size // 4
col_right = 3 * (size // 4)


def update(time):

    global h
    global step
    h_prev = h
    ex_prev = ex
    ey_prev = ey

    # update h
    h = h_prev + mu_arr * ((ex_prev[1:] - ex_prev[:-1]) - (ey_prev[:, 1:] - ey_prev[:, :-1]))

    # override an h column with a plane wave
    if 0 < time < pulseMid * 5:
        magnitude = (-time + pulseMid) * (1 / (sigma_t * math.sqrt(2 * math.pi))) * (
            math.exp(-((time - pulseMid) ** 2) / (2 * (sigma_t ** 2))))
        h[:, size//3] = magnitude
        h[:, 1 + (size//3)] += mu_arr[:, 1 + (size // 3)] * magnitude

    # update ex and ey, notice that we do not update the top or bottom row, first or last column
    ex[1:-1, : ] = ex_prev[1:-1, :] + eps_arr[1:, : ] * (h[1:, : ] - h[:-1, : ])
    ey[ :, 1:-1] = ey_prev[:, 1:-1] - eps_arr[ :, 1:] * (h[ :, 1:] - h[ :, :-1])

    if 0 < time < pulseMid * 5:
        magnitude = (-time + pulseMid) * (1 / (sigma_t * math.sqrt(2 * math.pi))) * (
            math.exp(-((time - pulseMid) ** 2) / (2 * (sigma_t ** 2))))

        # go left
        # ey[ :, 1 + (size // 3)] -= eps_arr[:, 1 + (size // 3)] * magnitude

        # go right
        ey[ :, size // 3] += eps_arr[ :, size // 3] * magnitude

    # absorption boundaries
    # h[0, :] = h[0, :] * (1 - (c * dt / ds)) + h[1, :] * (c * dt / ds)
    ex[0, :] = ex_prev[0, :] * (1 - (c * dt / ds)) + ex_prev[1, :] * (c * dt / ds)
    ey[0, : ] = ey_prev[0, : ] * (1 - (c * dt / ds)) + ey_prev[1, : ] * (c * dt / ds)
    # h[-1, :] = h[-1, :] * (1 - (c * dt / ds)) + h[-2, :] * (c * dt / ds)
    ex[-1, :] = ex_prev[-1, :] * (1 - (c * dt / ds)) + ex_prev[-2, :] * (c * dt / ds)
    ey[-1, : ] = ey_prev[-1, : ] * (1 - (c * dt / ds)) + ey_prev[-2, : ] * (c * dt / ds)

    step += 1

    return h


# plot properties
im = plt.imshow(h, animated=True, vmax=maxima, vmin=-maxima, aspect='auto')
plt.set_cmap('seismic')
ax.xaxis.set_ticks_position('top') # the rest is the same
ax.xaxis.set_label_position('top')
fig.colorbar(im)


def animate(time):
    if step % 10 == 0:
        ax.set_title("Time Step = {}".format(step))
        # plt.savefig('figs/' + str(int(step/50)) + '.png')

    im.set_array(update(time))
    return im,


anim = animation.FuncAnimation(fig, animate, frames=np.arange(0, end_time, dt), interval=1, blit=False, repeat=False)
# anim.save(filename="movie.gif", fps=20)

# Set up formatting for the movie files
# Writer = animation.writers['pillow']
# writer = Writer(metadata=dict(artist='Me'), bitrate=1800)
# anim.save('movie.gif', writer=writer)
plt.show()
exit()