import time
from analyser import FileLoader
from solver import Pulse
import matplotlib.pyplot as plt
import numpy as np

fileLoad = FileLoader("d_slit_test")

data_collectors = []

got_frequencies = False
frequency_axis = None
input_pulse = Pulse(fileLoad.pulses[0]['sigma_w'], fileLoad.dt, 0, fileLoad.pulses[0]['type'], fileLoad.end_time, omega_0=fileLoad.pulses[0]['omega_0'], step_frequency=fileLoad.step_frequency)

# input_pulse.plot_frequency_fft()
input_data = fileLoad.create_data_collector((2, 0.2))

for i in range(len(fileLoad.get_matrix()[0]) // 2):
    vertical_pos = 2 * i * fileLoad.ds
    data = fileLoad.create_data_collector((vertical_pos, 3.9))
    data.fft()

    if not got_frequencies:
        frequency_axis = data.fft_frequencies
        got_frequencies = True
    data.plot_frequency(show=False)
    data_collectors.append(data.fft_amplitude)
print(len(frequency_axis))
input_pulse.plot_frequency_fft()

plt.figure()
plt.title("Double slit diffraction")
plt.ylabel("Position (metres)")
plt.xlabel("Frequency (rad/s)")
im = plt.imshow(data_collectors, extent=[np.min(frequency_axis), np.max(frequency_axis), 4, 0], aspect='auto')
plt.show()
