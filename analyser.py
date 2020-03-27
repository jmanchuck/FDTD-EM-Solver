import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import json


class FileLoader:

    def __init__(self, file_name):
        self.file_name = file_name
        self.constants = None
        self.h_matrix = None
        self.size = None
        self.load()

    def load(self):
        print("Loading npy file...")
        self.h_matrix = np.load(self.file_name + ".npy")

        try:
            with open(self.file_name + ".txt") as json_file:
                print("Loading txt file...")
                self.constants = json.load(json_file)
                self.unpack_constants_from_json()
        except FileNotFoundError:
            print("Simulation constants txt file could not be found")

    def unpack_constants_from_json(self):
        self.dt = self.constants['dt']
        self.ds = self.constants['ds']
        self.length_x = self.constants['length_x']
        self.length_y = self.constants['length_y']
        self.stability = self.constants['stability']
        self.omega_max = self.constants['omega_max']

        self.end_time = self.constants['end_time']
        self.step_frequency = self.constants['step_frequency']

        self.material = np.array(self.constants['material'])
        self.size = len(self.material)

        self.pulses = self.constants['pulses']

    def play(self, interval=100, colourdepth=0.8):
        darkness_factor = 1 - colourdepth
        fig, ax = plt.subplots()

        def animate(i):
            matrix.set_array(self.h_matrix[i])
            fig.suptitle("Time step: {}".format(i * self.step_frequency))

        matrix = ax.imshow(self.h_matrix[0], vmax=darkness_factor*np.max(self.h_matrix), vmin=-darkness_factor*np.max(self.h_matrix), extent=[0, self.length_x, self.length_y, 0], cmap='seismic')
        plt.colorbar(matrix)
        plt.ylabel("Y (m)")
        plt.xlabel("X (m)")
        ani = animation.FuncAnimation(fig, animate, frames=len(self.h_matrix), interval=interval, repeat=False)
        plt.show()

    def plot_matrix_at_time(self, time):
        plt.figure()

        time_step = int(time / (self.dt * self.step_frequency))

        im = plt.imshow(self.h_matrix[time_step], vmax = np.max(self.h_matrix), vmin=-np.max(self.h_matrix),  extent=[0, self.length_x, self.length_y, 0], cmap='seismic')
        plt.colorbar(im)
        plt.suptitle("Time: " + str(time))
        plt.ylabel("Y (m)")
        plt.xlabel("X (m)")
        plt.show()

    def get_matrix(self):
        return self.h_matrix

    def convert(self, value):
        return int(value / self.ds)

    def create_data_collector(self, location):
        data_collector = DataCollector(self.h_matrix, location, self.constants)
        return data_collector


class DataCollector:

    def __init__(self, matrix, location, constants):
        self.matrix = matrix

        self.location = location
        self.constants = constants

        self.data = []

        # fast fourier transform plots
        self.fft_amplitude = None          # the y axis of the plot
        self.fft_frequencies = None  # the x axis of the plot

        self.unpack_constants()
        self.collect_all()

    def unpack_constants(self):
        self.dt = self.constants['dt']
        self.ds = self.constants['ds']
        self.length_x = self.constants['length_x']
        self.length_y = self.constants['length_y']
        self.stability = self.constants['stability']
        self.omega_max = self.constants['omega_max']

        self.end_time = self.constants['end_time']
        self.step_frequency = self.constants['step_frequency']

        self.i_pos = self.convert(self.location[0])
        self.j_pos = self.convert(self.location[1])

        self.time_vector = np.arange(0, len(self.matrix)) * self.dt * self.step_frequency

    def convert(self, value):
        return int(value / self.ds)

    def collect_all(self):
        self.data = self.matrix[:, self.i_pos, self.j_pos]

    def plot_time(self, show=True):
        plt.suptitle("Detected pulse")
        plt.xlabel("Time (seconds)")
        plt.ylabel("Magnitude")
        plt.plot(self.time_vector, self.data)
        if show:
            plt.show()

    def fft(self):
        # The sample spacing is NOT every timestep, since data collection may have happened every few timesteps
        # Multiply by 2π for frequency in units rad/s
        padded_zero_vector = np.zeros(2 ** 14 - len(self.data))
        padded_magnitude_vector = np.concatenate((self.data, padded_zero_vector))

        self.fft_frequencies = 2 * np.pi * np.fft.fftfreq(len(padded_magnitude_vector), d=self.dt * self.step_frequency)

        # Frequency mask to plot frequencies larger than 0 and smaller than bandwidth
        freq_mask = np.logical_and(self.fft_frequencies > 0, self.fft_frequencies < self.omega_max)
        self.fft_amplitude = np.absolute(np.fft.fft(padded_magnitude_vector)[freq_mask] * self.dt * self.step_frequency) ** 2
        self.fft_frequencies = self.fft_frequencies[freq_mask]

    def plot_frequency(self, show=True):
        if self.fft_amplitude is None:
            self.fft()

        plt.plot(self.fft_frequencies, self.fft_amplitude)
        if show:
            plt.xlabel("Frequency (rad/s)")
            plt.show()
