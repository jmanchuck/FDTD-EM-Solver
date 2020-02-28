import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from analyser import DataCollector
import time

# Physical constants
C = 2.99 * (10 ** 8)
MU = 1.26 * (10 ** (-6))
EPSILON = 8.85 * (10 ** (-12))


class Solver:

    def __init__(self, points_per_wavelength, stability, eps_r_max, mu_r_max, simulation_size, simulation_time=1):
        self.points_per_wavelength = points_per_wavelength  # mesh points per wavelength
        self.stability = stability  # time mesh stability factor

        # Material dependent, default is vacuum
        self.n_max = math.sqrt(eps_r_max * mu_r_max)

        # Pulse dependent
        self.omega_max = 0
        self.lambda_min = None
        self.length_x, self.length_y = simulation_size, simulation_size  # size of simulation is set
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

        self.h_list = []

        # material matrix
        self.material = None

        self.eps_coefficient_arr = None
        self.mu_coefficient_arr = None

        # simulation
        self.reflect = False
        self.reflectors = []
        self.boundaries = [False, False, False, False]  # up, down, left, right: True=reflect, False=absorb
        self.pulses = []
        self.pulses_max = None

        # save file
        self.save_file_name = None
        self.save_file = False

        # list of data collector objects
        self.data_collectors = []

    def update_constants(self):
        # Derived constants, calculate from user input but stays constant throughout simulation
        self.lambda_min = math.pi * 2 * C / (self.n_max * self.omega_max)
        self.ds = self.lambda_min / self.points_per_wavelength
        self.dt = min(self.ds * self.stability / C, math.pi / self.omega_max)  # mesh stability and Nyquist criterion
        self.size = int(self.length_x / self.ds)

        print("Smallest wavelength: {}\nMatrix size: {}".format(self.lambda_min, self.ds, self.size))

        # Resize matrices
        self.h = np.zeros((self.size, self.size))
        self.ex = np.zeros((self.size + 1, self.size))
        self.ey = np.zeros((self.size, self.size + 1))

        self.steps = int(self.end_time / self.dt)

        print("Number of time steps:", self.steps)

        # material properties matrix
        self.eps_arr = (np.ones((self.size, self.size)) * EPSILON)
        self.mu_arr = (np.ones((self.size, self.size)) * MU)

    def convert(self, si_unit):
        return int((si_unit / self.length_x) * self.size)

    def create_material(self):
        self.material = MaterialArray(self.size, self.length_x, self.length_y)
        return self.material

    def set_reflect_square(self, upper_left, lower_right):
        upper_left_converted = [self.convert(i) for i in upper_left]
        lower_right_converted = [self.convert(j) for j in lower_right]
        self.reflect = True
        reflect_arr = [upper_left_converted[0], lower_right_converted[0], upper_left_converted[1],
                       lower_right_converted[1]]
        self.reflectors.append(reflect_arr)

    # manually set the boundaries to reflect, on default call it sets all boundaries to be reflective
    def set_reflect_boundaries(self, up=True, down=True, left=True, right=True):
        self.boundaries = [up, down, left, right]

    def add_oscillating_pulse(self, sigma_w, location, omega_0, start_time=0, direction=None):
        if 3 * sigma_w + omega_0 > self.omega_max:
            self.omega_max = 3 * sigma_w + omega_0
            self.update_constants()
        location = [self.convert(i) for i in location]
        new_pulse = Pulse(sigma_w=sigma_w, dt=self.dt, location=location, omega_0=omega_0, start_time=start_time,
                          type="oscillate", direction=direction)
        self.pulses.append(new_pulse)

    def add_gaussian_pulse(self, sigma_w, location, start_time=0, direction=None):
        if 3 * sigma_w > self.omega_max:
            self.omega_max = 3 * sigma_w
            self.update_constants()
        location = [self.convert(i) for i in location]
        new_pulse = Pulse(sigma_w=sigma_w, dt=self.dt, location=location, start_time=start_time, type="gd",
                          direction=direction)
        self.pulses.append(new_pulse)

    def ready(self):
        return len(self.pulses) > 0

    def update_reflect(self):
        for square in self.reflectors:
            # reflecting boundary for square h
            self.h[square[0]:square[1], square[2]:square[3]] = 0

    def update(self, time):

        h_prev = self.h
        ex_prev = self.ex
        ey_prev = self.ey

        # override h field for pulse
        for pulse in self.pulses:
            if pulse.get_start_time() < time < pulse.get_end_time():
                magnitude = pulse.magnitude(time)
                if pulse.get_direction() is None:
                    self.h[pulse.get_row()][pulse.get_col()] = magnitude

        self.h = h_prev + self.mu_coefficient_arr * ((ex_prev[1:] - ex_prev[:-1]) - (ey_prev[:, 1:] - ey_prev[:, :-1]))
        # # override update equation at the TF/SF boundary
        for pulse in self.pulses:
            if pulse.get_start_time() < time < pulse.get_end_time() and not pulse.get_direction() is None:
                if pulse.get_direction() == "right":
                    self.h[:, pulse.get_col()] = h_prev[:, pulse.get_col()] + self.mu_coefficient_arr[:, pulse.get_col()] * \
                                                 ((ex_prev[1:, pulse.get_col()] - ex_prev[:-1, pulse.get_col()]) -
                                                  (ey_prev[:, 1 + pulse.get_col()] -
                                                   (ey_prev[:, pulse.get_col()] + pulse.magnitude(time) * math.sqrt(
                                                   MU / EPSILON))))

                elif pulse.get_direction() == "left":
                    self.h[:, pulse.get_col()] += self.mu_coefficient_arr[:, pulse.get_col()] * magnitude
                elif pulse.get_direction() == "up":
                    self.h[1 + pulse.get_row(), :] -= self.mu_coefficient_arr[1 + pulse.get_row(), :] * magnitude
                elif pulse.get_direction() == "down":
                    self.h[pulse.get_row(), :] -= self.mu_coefficient_arr[pulse.get_row(), :] * magnitude

        # if we have reflecting squares
        if self.reflect:
            self.update_reflect()

        # update ex and ey, notice that we do not update the top and bottom row, first and last column
        self.ex[1:-1, :] = ex_prev[1:-1, :] + self.eps_coefficient_arr[1:, :] * (self.h[1:, :] - self.h[:-1, :])
        self.ey[:, 1:-1] = ey_prev[:, 1:-1] - self.eps_coefficient_arr[:, 1:] * (self.h[:, 1:] - self.h[:, :-1])

        # override update equation at the TF/SF boundary
        for pulse in self.pulses:
            if pulse.get_start_time() < time < pulse.get_end_time() and not pulse.get_direction() is None:
                if pulse.get_direction() == "right":
                    self.ey[:, pulse.get_col()] = ey_prev[:, pulse.get_col()] - self.eps_coefficient_arr[:, pulse.get_col()] * (
                            self.h[:, pulse.get_col()] - pulse.magnitude(time) - (self.h[:, pulse.get_col() - 1]))

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

        for collector in self.data_collectors:
            collector.collect(self.h)

        self.step += 1

        return self.h

    def assign_pulse_max(self):

        # get the maximum of all pulses
        pulse_max = self.pulses[0].get_maximum()
        for pulse in self.pulses:
            if pulse.get_maximum() > pulse_max:
                pulse_max = pulse.get_maximum()

        self.pulses_max = pulse_max

        print('Solver max:', pulse_max)

    def solve(self, realtime=True):

        if not self.ready():
            raise Exception('No pulse added')

        self.eps_coefficient_arr = (self.dt / self.ds) * (1 / self.material.eps_arr)
        self.mu_coefficient_arr = (self.dt / self.ds) * (1 / self.material.mu_arr)

        self.assign_pulse_max()

        # vector of time values
        frames = np.arange(0, self.end_time, self.dt)

        if realtime:
            self.plot_and_show(frames)

        if not realtime:
            self.plot_and_save(frames)

        if self.save_file:
            np.save(self.save_file_name, self.h_list)

    def plot_and_show(self, frames):
        # instantiate animation plotting variables
        fig, ax1 = plt.subplots(figsize=(6, 6))
        im = ax1.imshow(self.h, animated=True, vmax=self.pulses_max, vmin=-self.pulses_max, aspect='auto', cmap='seismic')
        ax1.set_aspect('equal', 'box')
        ax1.xaxis.set_ticks_position('top')
        ax1.xaxis.set_label_position('top')
        fig.colorbar(im, ax=ax1)

        def animate(time):
            if self.step % 100 == 0:
                fig.suptitle("Time Step = {}".format(self.step))
                if self.save_file:
                    self.h_list.append(self.h)
            im.set_array(self.update(time))
            return im,

        anim = animation.FuncAnimation(fig, animate, frames=frames, interval=1, blit=False, repeat=False)

        plt.show()

    def plot_and_save(self, frames):
        for time in frames:
            self.update(time)
            if self.step % 20 == 0 and self.save_file:
                self.h_list.append(self.h)

    def save(self, file_name):
        self.save_file_name = file_name
        self.save_file = True


class MaterialArray:

    def __init__(self, size, real_size_x, real_size_y):
        # change eps matrix and mu matrix into the form of the constant we use in calculation

        self.eps_arr = np.ones((size, size)) * EPSILON
        self.mu_arr = np.ones((size, size)) * MU
        self.size = size
        self.length_x = real_size_x
        self.length_y = real_size_y

    def convert(self, si_unit):
        return int((si_unit / self.length_x) * self.size)

    def set_material_rect(self, upper_left, lower_right, epsilon_rel, mu_rel=1, convert=True):
        # sets a rectangle when given upper left and lower right coordinates (inclusive of lower right)
        if convert:
            upper_left = [self.convert(i) for i in upper_left]
            lower_right = [self.convert(j) for j in lower_right]

        self.eps_arr[upper_left[0]:lower_right[0], upper_left[1]:lower_right[1] + 1] = epsilon_rel * EPSILON
        self.mu_arr[upper_left[0]:lower_right[0], upper_left[1]:lower_right[1] + 1] = mu_rel * MU

    def set_material_convex(self, center: tuple, radius, thickness, epsilon_rel, mu_rel=1, convert=True):
        """
        Special note:
        Increase radius to have a more curved lens, but choose thickness accordingly otherwise it will go
        out of range.

        In general: larger radius (curvature) requires smaller thickness
        """
        if convert:
            displacer = self.convert(radius - (thickness / 2))
            center_i, center_j = self.convert(center[0]), self.convert(center[1])
            radius = self.convert(radius)
        else:
            displacer = radius - thickness // 2
            center_i, center_j = center[0], center[1]

        for j in range(-radius, -displacer):
            for i in range(-radius, radius):
                if i ** 2 + j ** 2 < radius ** 2:
                    self.eps_arr[center_i + i][center_j + j + displacer] = epsilon_rel * EPSILON
                    self.mu_arr[center_i + i][center_j + j + displacer] = mu_rel * MU

        for j in range(displacer, radius):
            for i in range(-radius, radius):
                if i ** 2 + j ** 2 < radius ** 2:
                    self.eps_arr[center_i + i][center_j + j - displacer] = epsilon_rel * EPSILON
                    self.mu_arr[center_i + i][center_j + j - displacer] = mu_rel * MU

    def set_fixed_length_waveguide(self, upper_left, waveguide_length, curved_portion_ratio, thickness, epsilon_rel,
                                   mu_rel=1):
        rect_length = (waveguide_length * (1 - curved_portion_ratio)) / 2
        radius = curved_portion_ratio * waveguide_length / (math.pi / 2)
        self.set_waveguide(upper_left, thickness, rect_length, radius + thickness / 2, epsilon_rel, mu_rel=mu_rel,
                     convert=True)

    def set_quarter_circle(self, center, radius, epsilon_rel, mu_rel=1, thickness=0, convert=True):
        center_i, center_j = center[0], center[1]
        if convert:
            center_i, center_j = self.convert(center_i), self.convert(center_j)
            radius = self.convert(radius)
            thickness = self.convert(thickness)
        for i in range(radius):
            for j in range(radius):
                if thickness == 0:
                    if i ** 2 + j ** 2 < radius ** 2:
                        self.eps_arr[center_i - i][center_j + j] = epsilon_rel * EPSILON
                        self.mu_arr[center_i - i][center_j + j] = mu_rel * MU
                else:
                    if radius ** 2 > i ** 2 + j ** 2 > (radius - thickness) ** 2:
                        self.eps_arr[center_i - i][center_j + j] = epsilon_rel * EPSILON
                        self.mu_arr[center_i - i][center_j + j] = mu_rel * MU

    def set_waveguide(self, start_point, thickness, rect_length, radius, epsilon_rel, mu_rel=1, convert=True):

        start_i = start_point[0]
        start_j = start_point[1]
        if convert:
            start_i = self.convert(start_i)
            start_j = self.convert(start_j)
            thickness = self.convert(thickness)
            rect_length = self.convert(rect_length)
            radius = self.convert(radius)

        self.set_material_rect((start_i, start_j), (start_i + thickness, start_j + rect_length), epsilon_rel, mu_rel=mu_rel, convert=False)

        quarter_circle_center = (start_i + radius, start_j + rect_length)
        self.set_quarter_circle(quarter_circle_center, radius, epsilon_rel, thickness=thickness, mu_rel=mu_rel, convert=False)

        down_rect_start = (quarter_circle_center[0], quarter_circle_center[1] + radius - thickness)
        self.set_material_rect(down_rect_start, (down_rect_start[0] + rect_length, down_rect_start[1] + thickness), epsilon_rel, mu_rel=mu_rel, convert=False)

    def plot(self):
        fig, ax = plt.subplots()
        image = ax.imshow(self.eps_arr, vmax=(np.max(self.eps_arr)), vmin=0, aspect='auto')
        ax.set_aspect('equal', 'box')
        ax.set_title('Material (epsilon relative)')

        plt.show()


class Pulse:

    def __init__(self, sigma_w, dt, location, start_time, type, omega_0=0, direction=None):
        # Gaussian derivative properties
        TYPES = ["gd", "oscillate"]
        DIRECTIONS = ["up", "down", "left", "right", None]
        self.sigma_w = sigma_w
        self.sigma_t = 1 / self.sigma_w
        self.omega_0 = omega_0
        self.dt = dt
        self.pulseMid = 3 * self.sigma_t  # the middle of the pulse (approximately 3 standard deviations)
        self._start_time = start_time
        self._end_time = self._start_time + self.pulseMid * 8
        self._location = location
        self._type = type
        self._maximum = self.calculate_max()
        self._omega_max = 3 * self.sigma_w + self.omega_0
        self._direction = direction
        self._end_step = int(1 + (self._end_time // self.dt))

        if self._type not in TYPES:
            raise Exception("Pulse type not recognised: {}".format(self._type))

        if self._direction not in DIRECTIONS:
            raise Exception("Pulse direction not recognised: {}".format(self.get_direction))

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
        biggest = 0
        while t < self.pulseMid * 3:
            if self._type == "gd":
                new = (-t + self.pulseMid) * (1 / (self.sigma_t * math.sqrt(2 * math.pi))) * (
                    math.exp(-((t - self.pulseMid) ** 2) / (2 * (self.sigma_t ** 2))))
            elif self._type == "oscillate":
                new = math.cos(t * self.omega_0) * (1 / (self.sigma_t * math.sqrt(2 * math.pi))) * (
                    math.exp(-((t - self.pulseMid) ** 2) / (2 * (self.sigma_t ** 2))))
            else:
                raise Exception("Invalid pulse")

            # if turning point found
            if new > biggest:
                biggest = new
            t += self.dt
        return biggest

    def plot(self):
        """
        Plot the oscillating pulse
        """
        data = []
        envelope = []
        t = 0
        while t < self._end_time:
            if self._type == "gd":
                value = (-t + self.pulseMid) * (1 / (self.sigma_t * math.sqrt(2 * math.pi))) * (
                    math.exp(-((t - self.pulseMid) ** 2) / (2 * (self.sigma_t ** 2))))
            elif self._type == "oscillate":
                value = math.cos(t * self.omega_0) * (1 / (self.sigma_t * math.sqrt(2 * math.pi))) * (
                    math.exp(-((t - self.pulseMid) ** 2) / (2 * (self.sigma_t ** 2))))
                envelope.append((1 / (self.sigma_t * math.sqrt(2 * math.pi))) * (
                    math.exp(-((t - self.pulseMid) ** 2) / (2 * (self.sigma_t ** 2)))))
            t += self.dt
            data.append(value)
        plt.plot(data)
        if self._type == "oscillate":
            plt.plot(envelope)
        plt.show()

    ### GET METHODS ###
    def get_maximum(self):
        return self._maximum

    def get_start_time(self):
        return self._start_time

    def get_direction(self):
        return self._direction

    def get_end_time(self):
        return self._end_time

    def get_end_step(self):
        return self._end_step

    def get_omega_max(self):
        return self._omega_max

    def get_row(self):
        return self._location[0]

    def get_col(self):
        return self._location[1]

    def get_type(self):
        return self._type


def test():
    # initiate variables
    sigma_w = 1 * 10 ** 9  # frequency bandwidth
    omega_0 = 5 * 10 ** 9  # central frequency

    print(2.99 * (10 ** 8) / omega_0)
    s = 10  # mesh points per wavelength
    stability = 0.2  # time mesh stability factor

    # initiate solver with user input variables
    solver = Solver(points_per_wavelength=s, stability=stability, eps_r_max=4, mu_r_max=1, simulation_size=3,
                    simulation_time=2.5 * 10 ** (-8))

    # add pulses (a pulse must be added to run the simulation)
    solver.add_oscillating_pulse(sigma_w, (0.8, 0.4), omega_0)  # this adds a point pulse

    mat = solver.create_material()

    mat.set_material_rect((2, 2), (2.8, 2.8), 3)

    mat.plot()

    exit()

    # solver.save('main_codetest')

    start_time = time.time()
    solver.solve(realtime=False)
    print("Time taken: ", time.time() - start_time)

if __name__ == '__main__':
    test()
    exit()
