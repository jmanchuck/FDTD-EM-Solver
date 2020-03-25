import time
from analyser import FileLoader
from solver import Pulse
import matplotlib.pyplot as plt
import numpy as np

# fileLoad = FileLoader("far_spaced_double_slit")
fileLoad = FileLoader("d_slit_test")

data_collectors = []

got_frequencies = False
frequency_axis = None
input_pulse = Pulse(fileLoad.pulses[0]['sigma_w'], fileLoad.dt, 0, fileLoad.pulses[0]['type'], fileLoad.end_time, omega_0=fileLoad.pulses[0]['omega_0'], step_frequency=fileLoad.step_frequency)


data_middle = fileLoad.create_data_collector((2, 3.95))
data_middle.collect_all()
data_middle.fft()
data_middle.plot_frequency(show=True)
print(len(data_middle.fft_amplitude))

input_pulse.plot_frequency_fft()
exit(0)

for i in range(len(fileLoad.get_matrix()[0]) // 2):
    vertical_pos = 2 * i * fileLoad.ds
    data = fileLoad.create_data_collector((vertical_pos, 3.9))
    data.collect_all()
    data.fft()

    if not got_frequencies:
        frequency_axis = data.fft_frequencies
        got_frequencies = True
    data.plot_frequency()
    data_collectors.append(data.fft_amplitude)

plt.figure()
plt.title("2 slit diffraction")

print(len(frequency_axis))

plt.ylabel("Position (metres)")
plt.xlabel("Frequency (rad/s)")

im = plt.imshow(data_collectors, extent=[np.min(frequency_axis), np.max(frequency_axis), 4, 0], aspect='auto')

plt.show()
