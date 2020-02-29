import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time


class FileLoader:

    def __init__(self, file_name):
        self.file_name = file_name
        self.h_matrix = None
        self.load()

    def load(self):
        self.h_matrix = np.load(self.file_name + ".npy")

    def play(self, interval=100, colourdepth=0.8):

        darkness_factor = 1 - colourdepth

        def animate(i):
            matrix.set_array(self.h_matrix[i])

        fig, ax = plt.subplots()
        matrix = ax.imshow(self.h_matrix[0], vmax=darkness_factor*np.max(self.h_matrix), vmin=-darkness_factor*np.max(self.h_matrix), cmap='seismic')
        plt.colorbar(matrix)
        ani = animation.FuncAnimation(fig, animate, frames=len(self.h_matrix), interval=interval, repeat=False)
        plt.show()

    def get_matrix(self):
        return self.h_matrix


class DataCollector:

    """
    This is the Data Collector.

    You can collect data from ANY of the 3 matrices (Ex, Ey or Hz) at 1 specific location (specific index).
    This object is accessed from the solver part of the code and automatically checks if you've added
    Data Collectors. You can add them using the method under solver: add_data_collector.

    Iterate through the list of Data Collectors in the solver to see plots of your data or fourier transforms of it.
    """

    def __init__(self, matrix, i_pos, j_pos):
        self.matrix = matrix
        self.i_pos = i_pos
        self.j_pos = j_pos
        self.data = []

        # fast fourier transform plots
        self.fft_amplitude = None          # the y axis of the plot
        self.fft_frequencies = None  # the x axis of the plot

    def collect_between(self, time_step_start, time_step_end):
        self.data = self.matrix[time_step_start:time_step_end, self.i_pos, self.j_pos]

    def collect_all(self):
        self.data = self.matrix[:, self.i_pos, self.j_pos]

    def plot_time(self):
        # fig, ax = plt.subplots()
        plt.plot(self.data)
        print(len(self.data))
        plt.show()

    def fft(self):
        fft_values = np.fft.fft(self.data)
        self.fft_amplitude = np.absolute(fft_values)
        self.fft_frequencies = np.fft.fftfreq(len(self.data))

    def plot_frequency(self):
        self.fft()
        plot_mask = self.fft_frequencies > 0
        plt.plot(self.fft_frequencies[plot_mask], self.fft_amplitude[plot_mask])
        plt.show()

    def get_data(self):
        return self.data


def test():

    loaded = FileLoader("double_slit")
    data = DataCollector(loaded.get_matrix(), 200, 200)
    data.collect_all()
    data.plot_time()
    data.plot_frequency()


if __name__ == "__main__":
    test()
