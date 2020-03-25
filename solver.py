import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from analyser import DataCollector, FileLoader
import time
import sys
import json

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
        self.dt = None  # discretized unit of space and time

        # simulation time in realtime and steps, length in real length
        self.steps = None
        self.end_time = simulation_time
        self.size = None
        self.time_vector = None

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
        self.save_dict_json = dict()
        self.step_frequency = None

    def save_to_json(self):
        self.save_dict_json['dt'] = self.dt
        self.save_dict_json['ds'] = self.ds
        self.save_dict_json['length_x'] = self.length_x
        self.save_dict_json['length_y'] = self.length_y
        self.save_dict_json['stability'] = self.stability
        self.save_dict_json['omega_max'] = self.omega_max

        self.save_dict_json['end_time'] = self.end_time
        self.save_dict_json['step_frequency'] = self.step_frequency

        self.save_dict_json['pulses'] = []
        for pulse in self.pulses:
            self.save_dict_json['pulses'].append({
                'sigma_w': pulse.sigma_w,
                'omega_0': pulse.omega_0,
                'omega_max': pulse.get_omega_max(),
                'row': pulse.get_row(),
                'col': pulse.get_col(),
                'type': pulse.get_type()
            })

        with open(self.save_file_name + '.txt', 'w') as outfile:
            json.dump(self.save_dict_json, outfile)

    def update_constants(self):
        # Derived constants, calculate from user input but stays constant throughout simulation
        self.lambda_min = math.pi * 2 * C / (self.n_max * self.omega_max)
        self.ds = self.lambda_min / self.points_per_wavelength
        self.dt = min(self.ds * self.stability / C, math.pi / self.omega_max)  # mesh stability and Nyquist criterion
        self.size = int(self.length_x / self.ds)

        # Resize matrices
        self.h = np.zeros((self.size, self.size))
        self.ex = np.zeros((self.size + 1, self.size))
        self.ey = np.zeros((self.size, self.size + 1))

        self.steps = int(self.end_time / self.dt)

        print("Matrix size: {}".format(self.size))
        print("Number of time steps:", self.steps)

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
                          type="oscillate", sim_end_time=self.end_time, direction=direction)
        self.pulses.append(new_pulse)

        return new_pulse

    def add_gaussian_pulse(self, sigma_w, location, start_time=0, direction=None):
        if 3 * sigma_w > self.omega_max:
            self.omega_max = 3 * sigma_w
            self.update_constants()
        location = [self.convert(i) for i in location]
        new_pulse = Pulse(sigma_w=sigma_w, dt=self.dt, location=location, start_time=start_time, type="gd", sim_end_time=self.end_time,
                          direction=direction)
        self.pulses.append(new_pulse)

        return new_pulse

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
                magnitude = pulse.magnitude(self.step)
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
                                                   (ey_prev[:, pulse.get_col()] + pulse.magnitude(self.step) * math.sqrt(
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
                            self.h[:, pulse.get_col()] - pulse.magnitude(self.step) - (self.h[:, pulse.get_col() - 1]))

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

    def assign_pulse_max(self):

        # get the maximum of all pulses
        pulse_max = self.pulses[0].pulse_maxima()
        for pulse in self.pulses:
            if pulse.pulse_maxima() > pulse_max:
                pulse_max = pulse.pulse_maxima()

        self.pulses_max = pulse_max

        print('Solver max:', pulse_max)

    def solve(self, realtime=True, step_frequency=5):

        if not self.ready():
            raise Exception('No pulse added')

        self.eps_coefficient_arr = (self.dt / self.ds) * (1 / self.material.eps_arr)
        self.mu_coefficient_arr = (self.dt / self.ds) * (1 / self.material.mu_arr)

        self.assign_pulse_max()

        # vector of time values
        self.time_vector = np.arange(0, self.end_time, self.dt)

        self.step_frequency = step_frequency

        if realtime and not self.save_file:
            self.plot_and_show()

        else:
            if not self.save_file:
                raise Exception("Must create a save file name for non-realtime simulation")
            else:
                self.plot_and_save(step_frequency=step_frequency)
                np.save(self.save_file_name, self.h_list)
                self.save_to_json()

    def plot_and_show(self):
        # instantiate animation plotting variables
        fig, ax1 = plt.subplots(figsize=(6, 6))
        im = ax1.imshow(self.h, animated=True, vmax=self.pulses_max, vmin=-self.pulses_max, aspect='auto', cmap='seismic')
        ax1.set_aspect('equal', 'box')
        ax1.xaxis.set_ticks_position('top')
        ax1.xaxis.set_label_position('top')
        fig.colorbar(im, ax=ax1)

        def animate(time):
            if self.step % 20 == 0:
                fig.suptitle("Time Step = {}".format(self.step))
                if self.save_file:
                    self.h_list.append(self.h)
            im.set_array(self.update(time))
            return im,

        anim = animation.FuncAnimation(fig, animate, frames=self.time_vector, interval=1, blit=False, repeat=False)

        plt.show()

    def plot_and_save(self, step_frequency):
        for time in self.time_vector:
            self.update(time)
            sys.stdout.write("\r%d%%" % (100 * self.step // self.steps))
            sys.stdout.flush()
            if self.step % step_frequency == 0 and self.save_file:
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

    def __init__(self, sigma_w, dt, start_time, type, sim_end_time, location=None, omega_0=0, direction=None, step_frequency=1):
        print("--- Pulse Info ---")
        TYPES = ["gd", "oscillate"]
        DIRECTIONS = ["up", "down", "left", "right", None]
        self.sigma_w = sigma_w
        self.sigma_t = 1 / self.sigma_w
        self.omega_0 = omega_0
        self.dt = dt
        self.pulseMid = 3 * self.sigma_t  # the middle of the pulse (approximately 3 standard deviations)

        self._start_time = start_time
        self._end_time = self._start_time + self.pulseMid * 3
        self._sim_end_time = sim_end_time
        self._location = location
        self._type = type
        self.step_frequency = step_frequency

        self._omega_max = 3 * self.sigma_w + self.omega_0
        self._direction = direction
        self._end_step = int(1 + (self._end_time // self.dt))

        print("Pulse type: {}, sigma_w: {}, omega_0: {}, omega_max: {}".format(self._type, self.sigma_w, self.omega_0, self._omega_max))

        self._time_vector = np.arange(0, self._sim_end_time, self.dt)
        self._magnitude_vector = None
        self.create_magnitude_vector()
        self._maximum = self.pulse_maxima()

        if self._type not in TYPES:
            raise Exception("Pulse type not recognised: {}".format(self._type))

        if self._direction not in DIRECTIONS:
            raise Exception("Pulse direction not recognised: {}".format(self.get_direction))

        if self._type == "oscillate" and not self.omega_0:
            raise Exception("Missing arguments f_0: {}".format(self.omega_0))

    def magnitude(self, timestep):
        return self._magnitude_vector[timestep]

    def pulse_maxima(self):
        return np.max(self._magnitude_vector)

    def create_magnitude_vector(self):
        t = self._time_vector
        if self._type == "gd":
            self._magnitude_vector = (-t + self.pulseMid) * (1 / (self.sigma_t * np.sqrt(2 * math.pi))) * (
                    np.exp(-((t - self.pulseMid) ** 2) / (2 * (self.sigma_t ** 2))))
        elif self._type == "oscillate":
            self._magnitude_vector = np.cos(t * self.omega_0) * (1 / (self.sigma_t * np.sqrt(2 * math.pi))) * (
                    np.exp(-((t - self.pulseMid) ** 2) / (2 * (self.sigma_t ** 2))))
        print("Pulse duration: {} seconds, {} steps".format(self._time_vector[-1], len(self._time_vector)))

    def plot_time(self, scaled_plot=True, show=True):
        plt.plot(self._time_vector, self._magnitude_vector)
        if self._type == "oscillate":
            t = np.arange(0, self._end_time, self.dt)
            envelope = (1 / (self.sigma_t * np.sqrt(2 * math.pi))) * (
                    np.exp(-((t - self.pulseMid) ** 2) / (2 * (self.sigma_t ** 2))))
            plt.plot(self._time_vector, envelope)

        plt.suptitle("{} pulse".format(self._type))
        plt.xlabel("Time (seconds)")
        plt.ylabel("Magnitude")
        if scaled_plot:
            plt.xlim(0, self._end_time)

        if show:
            plt.show()

    def plot_frequency(self, show=True):
        plt.xlabel("Frequency")
        plt.suptitle("{} pulse".format(self._type))
        frequencies_vector = np.linspace(self._omega_max, 1000)
        if self._type == "oscillate":
            plot_vector = self.sigma_t * 0.5 * (np.exp((-(frequencies_vector - self.omega_0) ** 2 / (2 * self.sigma_w ** 2))) +
                                                np.exp((-(frequencies_vector + self.omega_0) ** 2 / (2 * self.sigma_w ** 2))))
            plt.plot(frequencies_vector, plot_vector)
        elif self._type == "gd":
            plot_vector = np.abs(frequencies_vector) * self.sigma_t * np.exp(-((frequencies_vector ** 2) / (2 * self.sigma_w ** 2)))
            plt.plot(frequencies_vector, plot_vector)

        if show:
            plt.show()

    def plot_frequency_fft(self, show=True):
        magnitude_vector = self._magnitude_vector[::self.step_frequency]

        # Provides number of samples and their spacing in time domain
        # Multiply by 2π to obtain frequency axis in units of rad/s
        padded_zero_vector = np.zeros(2 ** 14 - len(magnitude_vector))
        print("Input pulse magnitude vector length:", len(magnitude_vector))
        padded_magnitude_vector = np.concatenate((magnitude_vector, padded_zero_vector))
        fft_frequencies = 2 * np.pi * np.fft.fftfreq(len(padded_magnitude_vector), d=self.dt * self.step_frequency)
        fft_values = np.fft.fft(padded_magnitude_vector)
        freq_mask = np.logical_and(fft_frequencies > 0, fft_frequencies < self._omega_max)
        print("Number of frequency points:", len(fft_frequencies[freq_mask]))

        # multiply by 2π for changing the units to rad / s instead of Hz
        plt.plot(fft_frequencies[freq_mask], np.abs(fft_values[freq_mask]) ** 2)

        if show:
            plt.suptitle("{} pulse (FFT)".format(self._type))
            plt.xlabel("Frequency (rad/s)")
            plt.show()

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
    stability = 0.1  # time mesh stability factor

    # initiate solver with user input variables
    solver = Solver(points_per_wavelength=s, stability=stability, eps_r_max=4, mu_r_max=1, simulation_size=3,
                    simulation_time=2.5 * 10 ** (-8))

    # add pulses (a pulse must be added to run the simulation)
    pulse = solver.add_oscillating_pulse(sigma_w, (0.8, 0.4), omega_0)  # this adds a point pulse

    pulse.plot_time()
    pulse.plot_frequency()
    pulse.plot_frequency_fft()

    mat = solver.create_material()

    mat.set_material_rect((2, 2), (2.8, 2.8), 3)

    solver.solve()


if __name__ == '__main__':
    test()
    exit()
