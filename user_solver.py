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

    def __init__(self, s, stability, simulation_time=1):
        self.s = s  # mesh points per wavelength
        self.stability = stability  # time mesh stability factor

        # Material dependent, default is vacuum
        self.n_max = math.sqrt(1 * 1)

        # Pulse dependent
        self.omega_max = 0
        self.lambda_min = None
        self.l_x, self.l_y = None, None  # size of simulation is 50 wavelengths
        self.ds = None
        self.dt = None  # mesh stability and Nyquist criterion

        # simulation time in realtime and steps, length in real length
        self.steps = None
        self.end_time = simulation_time
        self.size = None

        # simulation variables
        self.time = 0
        self.step = 0
        self.h = None
        self.ex = None
        self.ey = None
        self.entropy = []

        # material properties matrix
        self.eps_arr = None
        self.mu_arr = None

        # simulation
        self.reflect = False
        self.reflectors = []
        self.boundaries = [False, False, False, False]  # up, down, left, right: True=reflect, False=absorb
        self.pulses = []
        self.pulses_max = 0

    def update_constants(self):
        # Derived constants, calculate from user input but stays constant throughout simulation

        self.lambda_min = math.pi * 2 * C / (self.n_max * self.omega_max)
        self.l_x, self.l_y = 30 * self.lambda_min, 30 * self.lambda_min  # size of simulation is 50 wavelengths
        self.ds = self.lambda_min / self.s
        self.dt = min(self.ds * self.stability / C, math.pi / self.omega_max)  # mesh stability and Nyquist criterion
        self.size = int(self.l_x / self.ds)

        # Resize matrices
        self.h = np.zeros((self.size, self.size))
        self.ex = np.zeros((self.size + 1, self.size))
        self.ey = np.zeros((self.size, self.size + 1))

        self.steps = self.end_time / self.dt

        # material properties matrix
        self.eps_arr = (self.dt / self.ds) * (np.ones((self.size, self.size)) / EPSILON)
        self.mu_arr = (self.dt / self.ds) * (np.ones((self.size, self.size)) / MU)

    def add_material_square(self, upper_left, lower_right, epsilon_rel=1, mu_rel=1):
        # adding square material in top right by default for now
        self.eps_arr[upper_left[0]:lower_right[0], upper_left[1]:lower_right[1]] *= epsilon_rel
        self.mu_arr[upper_left[0]:lower_right[0], upper_left[1]:lower_right[1]] *= mu_rel

        self.n_max = math.sqrt((np.max(self.eps_arr) / EPSILON) * (np.max(self.mu_arr) / MU))

    def add_reflect_square(self, upper_left, lower_right):
        self.reflect = True
        reflect_arr = [upper_left[0], lower_right[0], upper_left[1], lower_right[1]]
        self.reflectors.append(reflect_arr)

    # manually set the boundaries to reflect, on default call it sets all boundaries to be reflective
    def set_reflect_boundaries(self, up=True, down=True, left=True, right=True):
        self.boundaries = [up, down, left, right]

    def add_oscillating_pulse(self, sigma_w, location, omega_0, start_time=0, direction=None):
        if 3 * sigma_w + omega_0 > self.omega_max:
            self.omega_max = 3 * sigma_w + omega_0
            self.update_constants()
        new_pulse = Pulse(sigma_w=sigma_w, dt=self.dt, location=location, omega_0=omega_0, start_time=start_time,
                          type="oscillate", direction=direction)
        self.load_pulse(new_pulse)

    def add_gaussian_pulse(self, sigma_w, location, start_time=0, direction=None):
        if 3 * sigma_w > self.omega_max:
            self.omega_max = 3 * sigma_w
            self.update_constants()
        new_pulse = Pulse(sigma_w=sigma_w, dt=self.dt, location=location, start_time=start_time, type="gd",
                          direction=direction)
        self.load_pulse(new_pulse)

    def load_pulse(self, pulse):
        self.pulses.append(pulse)
        if pulse.maximum() > self.pulses_max:
            self.pulses_max = pulse.maximum()

    def ready(self):
        return len(self.pulses) > 0

    def update(self, time):

        h_prev = self.h
        ex_prev = self.ex
        ey_prev = self.ey

        # override h field for pulse
        for pulse in self.pulses:
            if pulse.start_time() < time < pulse.end_time():
                magnitude = pulse.magnitude(time)
                if pulse.direction() is None:
                    self.h[pulse.row()][pulse.col()] = magnitude

        self.h = h_prev + self.mu_arr * ((ex_prev[1:] - ex_prev[:-1]) - (ey_prev[:, 1:] - ey_prev[:, :-1]))

        # override update equation at the TF/SF boundary
        for pulse in self.pulses:
            if pulse.start_time() < time < pulse.end_time() and not pulse.direction() is None:
                if pulse.direction() == "right":
                    # self.h[:, pulse.col()] = h_prev[:, pulse.col()] + self.mu_arr[:, pulse.col()] * \
                    #                          ((ex_prev[1:, pulse.col()] - ex_prev[:-1, pulse.col()]) -
                    #                           (ey_prev[:, 1 + pulse.col()] -
                    #                            (ey_prev[:, pulse.col()] + pulse.magnitude(time) * math.sqrt(MU / EPSILON))))
                    self.h[:, pulse.col()] += self.mu_arr[:, pulse.col()] * pulse.magnitude(time) * math.sqrt(MU / EPSILON)

                elif pulse.direction() == "left":
                    self.h[:, pulse.col()] += self.mu_arr[:, pulse.col()] * magnitude
                elif pulse.direction() == "up":
                    self.h[1 + pulse.row(), :] -= self.mu_arr[1 + pulse.row(), :] * magnitude
                elif pulse.direction() == "down":
                    self.h[pulse.row(), :] -= self.mu_arr[pulse.row(), :] * magnitude

        # if we have reflecting squares
        if self.reflect:
            for square in self.reflectors:
                # reflecting boundary for square h
                self.h[square[0]:square[1], square[2]:square[3]] = 0

        # update ex and ey, notice that we do not update the top and bottom row, first and last column
        self.ex[1:-1, :] = ex_prev[1:-1, :] + self.eps_arr[1:, :] * (self.h[1:, :] - self.h[:-1, :])
        self.ey[:, 1:-1] = ey_prev[:, 1:-1] - self.eps_arr[:, 1:] * (self.h[:, 1:] - self.h[:, :-1])

        # override update equation at the TF/SF boundary
        for pulse in self.pulses:
            if pulse.start_time() < time < pulse.end_time() and not pulse.direction() is None:
                if pulse.direction() == "right":
                    self.ey[:, pulse.col()] = ey_prev[:, pulse.col()] - self.eps_arr[:, pulse.col()] * (
                                self.h[:, pulse.col()] - pulse.magnitude(time) - (self.h[:, pulse.col() - 1]))
                elif pulse.direction() == "left":
                    self.ey[:, pulse.col()] += self.eps_arr[:, pulse.col()] * self.h[:, pulse.col()]
                elif pulse.direction() == "up":
                    self.ex[1 + pulse.row(), :] -= self.eps_arr[1 + pulse.row(), :] * magnitude
                elif pulse.direction() == "down":
                    self.ex[pulse.row(), :] += self.eps_arr[pulse.row(), :] * magnitude

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

        # overriding values for reflecting squares or objects that were added
        if self.reflect:
            for square in self.reflectors:
                # reflecting boundary for square ex and ey
                self.ex[square[0]:square[1], square[2]:square[3]] = 0
                self.ey[square[0]:square[1], square[2]:square[3]] = 0

        self.step += 1

        return self.h

    def solve(self):

        if not self.ready():
            raise Exception('No pulse added')
        # instantiate animation plotting variables
        frames = np.arange(0, self.end_time, self.dt)
        fig, ax1 = plt.subplots(figsize=(6, 6))

        im = ax1.imshow(self.h, animated=True, vmax=self.pulses_max, vmin=-self.pulses_max, aspect='auto',
                        cmap='seismic')
        ax1.set_aspect('equal', 'box')
        ax1.xaxis.set_ticks_position('top')  # the rest is the same
        ax1.xaxis.set_label_position('top')
        fig.colorbar(im, ax=ax1)

        def animate(time):
            if self.step % 10 == 0:
                fig.suptitle("Time Step = {}".format(self.step))
                # plt.savefig('figs/' + str(int(step/50)) + '.png')

            im.set_array(self.update(time))
            return im,

        anim = animation.FuncAnimation(fig, animate, frames=frames, interval=1, blit=False, repeat=False)

        plt.show()


# this class should be hidden from the user, user should NOT be able to access this class
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

    def __init__(self, sigma_w, dt, location, start_time, type, omega_0=0, direction=None):
        # Gaussian derivative properties
        TYPES = ["gd", "oscillate"]
        DIRECTIONS = ["up", "down", "left", "right", None]
        self.sigma_w = sigma_w
        self.sigma_t = 1 / self.sigma_w
        self.dt = dt
        self.pulseMid = 3 * self.sigma_t  # the middle of the gaussian derivative (approximate)
        self._start_time = start_time
        self._end_time = self._start_time + self.pulseMid * 8
        self._location = location
        self._type = type
        self._maximum = self.calculate_max()
        self.omega_0 = omega_0
        self._omega_max = 3 * self.sigma_w + self.omega_0
        self._direction = direction
        self._end_step = int(1 + (self._end_time // self.dt))

        if self._type not in TYPES:
            raise Exception("Pulse type not recognised: {}".format(self._type))

        if self._direction not in DIRECTIONS:
            raise Exception("Pulse direction not recognised: {}".format(self.direction))

        if self._type == "oscillate" and not self.omega_0:
            raise Exception("Missing arguments f_0: {}".format(self.omega_0))

    def magnitude(self, time):
        # shift the calculation
        time -= self._start_time
        if self._type == "gd":
            return (-time + self.pulseMid) * (1 / (self.sigma_t * math.sqrt(2 * math.pi))) * (
                math.exp(-((time - self.pulseMid) ** 2) / (2 * (self.sigma_t ** 2))))
        if self._type == "oscillate":
            return math.cos(time * self.omega_0) * (1 / (self.sigma_t * math.sqrt(2 * math.pi))) * (
                math.exp(-((time - self.pulseMid) ** 2) / (2 * (self.sigma_t ** 2))))

    def calculate_max(self):
        """
        :return: the maximum of the pulse via calculating values until turning point
        """
        t = 0
        current = 0
        prev = 0
        while t < 2 * self.pulseMid:
            if self._type == "gd":
                current = (-t + self.pulseMid) * (1 / (self.sigma_t * math.sqrt(2 * math.pi))) * (
                    math.exp(-((t - self.pulseMid) ** 2) / (2 * (self.sigma_t ** 2))))
            elif self._type == "oscillate":
                current = (1 / (self.sigma_t * math.sqrt(2 * math.pi))) * (
                    math.exp(-((t - self.pulseMid) ** 2) / (2 * (self.sigma_t ** 2))))

            # if turning point found
            if current < prev:
                return prev
            t += self.dt
            prev = current

        raise Exception("Could not find maximum")

    def maximum(self):
        return self._maximum

    def start_time(self):
        return self._start_time

    def direction(self):
        return self._direction

    def end_time(self):
        return self._end_time

    def end_step(self):
        return self._end_step

    def omega_max(self):
        return self._omega_max

    def row(self):
        return self._location[0]

    def col(self):
        return self._location[1]

    def type(self):
        return self._type

    def plot(self):
        data = []
        t = 0
        while t < 5 * self.pulseMid:
            if self._type == "gd":
                value = (-t + self.pulseMid) * (1 / (self.sigma_t * math.sqrt(2 * math.pi))) * (
                    math.exp(-((t - self.pulseMid) ** 2) / (2 * (self.sigma_t ** 2))))
            elif self._type == "oscillate":
                value = math.cos(t * self.omega_0) * (1 / (self.sigma_t * math.sqrt(2 * math.pi))) * (
                    math.exp(-((t - self.pulseMid) ** 2) / (2 * (self.sigma_t ** 2))))
            t += self.dt
            data.append(value)
        plt.plot(data)
        plt.show()
        exit()


# improvements:
# - take into account epsilon max to determine how we will sample a wavelength (to get dx or ds)
# - user inputs real time simulation length in seconds, use dt to calculate how many steps


def test():
    # initiate variables
    sigma_w = 1  # frequency bandwidth
    omega_0 = 3  # central frequency
    s = 10  # mesh points per wavelength
    stability = 0.2  # time mesh stability factor

    # initiate solver with user input variables
    solver = Solver(s=s, stability=stability, simulation_time=1000)

    # add pulse
    # solver.add_oscillating_pulse(sigma_w, (200, 150), omega_0, direction="right")

    solver.add_gaussian_pulse(sigma_w, (200, 90), direction="right")

    # add a reflecting square in top left
    solver.add_reflect_square((3 * solver.size // 8, 3 * solver.size // 8), (5 * solver.size // 8, 5 * solver.size // 8))

    # change the boundary to have reflect on the bottom
    solver.set_reflect_boundaries(up=True, down=True, left=False, right=True)

    print(solver.steps)
    print(solver.pulses[0].end_step())
    # exit()

    # add a material in the bottom right
    # solver.add_material_square((3 * solver.size // 4, 3 * solver.size // 4), (4 * solver.size // 5, 4 * solver.size // 5),
    #                            epsilon_rel=8)

    solver.solve()


if __name__ == '__main__':
    test()
    exit()

    # EXAMPLE INPUTS
    sigma_w = 1  # frequency bandwidth
    omega_0 = 0  # central frequency
    s = 10  # mesh points per wavelength
    stability = 0.2  # time mesh stability factor
    solver = Solver(s=s, stability=stability, simulation_time=1000)
    solver.set_reflect_boundaries(down=False, right=False)
    solver.add_gaussian_pulse(sigma_w, (50, 50))
    solver.add_oscillating_pulse(sigma_w, (180, 180), omega_0=3, start_time=10)
    solver.add_material_square((2 * solver.size // 3, 0), (-1, solver.size // 3), epsilon_rel=5, mu_rel=0.8)
    solver.solve()
    exit()
