import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from time import process_time

# Physical constants
C = 2.99 * (10 ** 8)
MU = 1.26 * (10 ** (-6))
EPSILON = 8.85 * (10 ** (-12))


class Solver:

    def __init__(self, sigma_w, omega_0, s, stability):
        self.sigma_w = sigma_w
        self.omega_0 = omega_0  # central frequency
        self.s = s  # mesh points per wavelength
        self.stability = stability  # time mesh stability factor

        self.epsilon_relative = 2  # change this to a list in the future
        self.mu_relative = 1  # change this to a list

        # Derived constants, calculate from user input but stays constant throughout simulation
        self.n_max = math.sqrt(self.epsilon_relative * self.mu_relative)
        self.sigma_t = 1 / self.sigma_w
        self.omega_max = self.omega_0 + 3 * self.sigma_w
        self.lambda_min = math.pi * 2 * C / (self.n_max * self.omega_max)
        self.l_x, self.l_y = 30 * self.lambda_min, 30 * self.lambda_min  # size of simulation is 50 wavelengths
        self.ds = self.lambda_min / s
        self.dt = min(self.ds * stability / C, (4 * math.pi) / self.omega_max)

        # simulation time in realtime and steps, length in real length
        self.steps = 1500
        self.end_time = self.dt * self.steps
        self.size = int(self.l_x / self.ds)

        # simulation variables
        self.time = 0
        self.step = 0
        self.h = np.zeros((self.size, self.size))
        self.ex = np.zeros((self.size + 1, self.size))
        self.ey = np.zeros((self.size, self.size + 1))

        # animation plotting
        self.frames = np.arange(0, self.end_time, self.dt)
        print(self.frames)
        self.fig, self.ax = plt.subplots()

        # material properties matrix
        self.eps_arr = (self.dt / self.ds) * (np.ones((self.size, self.size)) / EPSILON)
        self.mu_arr = (self.dt / self.ds) * (np.ones((self.size, self.size)) / MU)

        # simulation
        # self.ready = False
        self.reflect = False
        self.reflectors = []
        self.boundaries = [False, False, False, False]  # up, down, left, right: True=reflect, False=absorb
        self.pulse = Pulse(self.sigma_t, self.dt)
        self.pulse_max = self.pulse.maxima()
        self.pulse_pos = [self.size // 2, self.size // 2]

    def add_material_square(self, upper_left, lower_right, epsilon_rel=1, mu_rel=1):
        # adding square material in top right by default for now
        self.eps_arr[upper_left[0]:lower_right[0], upper_left[1]:lower_right[1]] *= epsilon_rel
        self.mu_arr[upper_left[0]:lower_right[0], upper_left[1]:lower_right[1]] *= mu_rel

    def add_reflect_square(self, upper_left, lower_right):
        self.reflect = True
        reflect_arr = [upper_left[0], lower_right[0], upper_left[1], lower_right[1]]
        self.reflectors.append(reflect_arr)

    # manually set the boundaries to reflect, on default call it sets all boundaries to be reflective
    def set_boundaries(self, up=True, down=True, left=True, right=True):
        self.boundaries = [up, down, left, right]

    def update(self, time):
        h_prev = self.h
        ex_prev = self.ex
        ey_prev = self.ey

        self.h = h_prev + self.mu_arr * ((ex_prev[1:] - ex_prev[:-1]) - (ey_prev[:, 1:] - ey_prev[:, :-1]))

        # update h pulse here

        # override h bottom left with Gaussian
        if 0 < time < self.pulse.end_time:
            magnitude = self.pulse.magnitude(time)
            self.h[self.pulse_pos[0]][self.pulse_pos[1]] = magnitude
            # print(magnitude)

        # if we have reflecting squares
        if self.reflect:
            for square in self.reflectors:
                # reflecting boundary for square h
                self.h[square[0]:square[1], square[2]:square[3]] = 0

        # update ex and ey, notice that we do not update the top and bottom row, first and last column
        self.ex[1:-1, :] = ex_prev[1:-1, :] + self.eps_arr[1:, :] * (self.h[1:, :] - self.h[:-1, :])
        self.ey[:, 1:-1] = ey_prev[:, 1:-1] - self.eps_arr[:, 1:] * (self.h[:, 1:] - self.h[:, :-1])

        # if the boundary is NOT reflective, apply absorb equation, order is up, down, left and right borders
        if not self.boundaries[0]:
            self.h[0, :] = h_prev[0, :] * (1 - (C * self.dt / self.ds)) + h_prev[1, :] * (C * self.dt / self.ds)
            self.ex[0, :] = ex_prev[0, :] * (1 - (C * self.dt / self.ds)) + ex_prev[1, :] * (C * self.dt / self.ds)
            self.ey[0, :] = ey_prev[0, :] * (1 - (C * self.dt / self.ds)) + ey_prev[1, :] * (C * self.dt / self.ds)
        if not self.boundaries[1]:
            self.h[-1, :] = h_prev[-1, :] * (1 - (C * self.dt / self.ds)) + h_prev[-2, :] * (C * self.dt / self.ds)
            self.ex[-1, :] = ex_prev[-1, :] * (1 - (C * self.dt / self.ds)) + ex_prev[-2, :] * (C * self.dt / self.ds)
            self.ey[-1, :] = ey_prev[-1, :] * (1 - (C * self.dt / self.ds)) + ey_prev[-2, :] * (C * self.dt / self.ds)
        if not self.boundaries[2]:
            self.h[:, 0] = h_prev[:, 0] * (1 - (C * self.dt / self.ds)) + h_prev[:, 1] * (C * self.dt / self.ds)
            self.ex[:, 0] = ex_prev[:, 0] * (1 - (C * self.dt / self.ds)) + ex_prev[:, 1] * (C * self.dt / self.ds)
            self.ey[:, 0] = ey_prev[:, 0] * (1 - (C * self.dt / self.ds)) + ey_prev[:, 1] * (C * self.dt / self.ds)
        if not self.boundaries[3]:
            self.h[:, -1] = h_prev[:, -1] * (1 - (C * self.dt / self.ds)) + h_prev[:, -2] * (C * self.dt / self.ds)
            self.ex[:, -1] = ex_prev[:, -1] * (1 - (C * self.dt / self.ds)) + ex_prev[:, -2] * (C * self.dt / self.ds)
            self.ey[:, -1] = ey_prev[:, -1] * (1 - (C * self.dt / self.ds)) + ey_prev[:, -2] * (C * self.dt / self.ds)

        if self.reflect:
            for square in self.reflectors:
                # reflecting boundary for square ex and ey
                self.ex[square[0]:square[1], square[2]:square[3]] = 0
                self.ey[square[0]:square[1], square[2]:square[3]] = 0

        self.step += 1

        return self.h

    def solve(self):

        def animate(time):
            if self.step % 10 == 0:
                self.ax.set_title("Time Step = {}".format(self.step))
                # plt.savefig('figs/' + str(int(step/50)) + '.png')

            im.set_array(self.update(time))
            return im,

        im = plt.imshow(self.h, animated=True, vmax=self.pulse_max, vmin=-self.pulse_max, aspect='auto')
        plt.set_cmap('seismic')
        self.ax.xaxis.set_ticks_position('top')  # the rest is the same
        self.ax.xaxis.set_label_position('top')
        self.fig.colorbar(im)

        anim = animation.FuncAnimation(self.fig, animate, frames=self.frames, interval=1, blit=False, repeat=False)
        plt.show()
        exit()


class Pulse:
    """
    This is only for Gaussian derivative, in the future this will be
    Abstract class API with methods:
        magnitude
        maxima
        plot
    Inheriting classes:
        gaussian derivative
        oscillating pulse
        other types of pulses
    """

    def __init__(self, sigma_t, dt, type="gd", f_0=None):
        # Gaussian derivative properties
        self.sigma_t = sigma_t
        self.dt = dt
        self.pulseMid = 3 * sigma_t  # the middle of the gaussian derivative (approximate)
        self._maxima = self.calc_maxima()
        self.end_time = self.pulseMid * 5

    def magnitude(self, time):
        return (-time + self.pulseMid) * (1 / (self.sigma_t * math.sqrt(2 * math.pi))) * (
            math.exp(-((time - self.pulseMid) ** 2) / (2 * (self.sigma_t ** 2))))

    def calc_maxima(self):
        t = 0
        current = 0
        prev = 0
        while t < self.pulseMid:
            current = (-t + self.pulseMid) * (1 / (self.sigma_t * math.sqrt(2 * math.pi))) * (
                math.exp(-((t - self.pulseMid) ** 2) / (2 * (self.sigma_t ** 2))))
            if current < prev:
                return prev
            t += self.dt
            prev = current
        return current

    def maxima(self):
        return self._maxima

    def end_time(self):
        return self.end_time

    def plot(self):
        data = []
        t = 0
        while t < 5 * self.pulseMid:
            gauss = (-t + self.pulseMid) * (1 / (self.sigma_t * math.sqrt(2 * math.pi))) * (
                math.exp(-((t - self.pulseMid) ** 2) / (2 * (self.sigma_t ** 2))))
            t += self.dt
            data.append(gauss)
        plt.plot(data)
        plt.show()
        exit()


# improvements:
# - take into account epsilon max to determine how we will sample a wavelength (to get dx or ds)
# - user inputs real time simulation length in seconds, use dt to calculate how many steps

if __name__ == '__main__':
    # EXAMPLE INPUTS
    sigma_w = 1  # frequency bandwidth
    omega_0 = 0  # central frequency
    s = 10  # mesh points per wavelength
    stability = 0.2  # time mesh stability factor
    solver = Solver(sigma_w=sigma_w, omega_0=omega_0, s=s, stability=stability)
    solver.solve()
    exit()