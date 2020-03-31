from analyser import FileLoader
from solver import Pulse
import matplotlib.pyplot as plt
import numpy as np

fileLoad = FileLoader("double_slit_short_pulse")

data_collectors = []

got_frequencies = False
frequency_axis = None
input_pulse = Pulse(fileLoad.pulses[0]['sigma_w'], fileLoad.dt, 0, fileLoad.pulses[0]['type'], fileLoad.end_time, omega_0=fileLoad.pulses[0]['omega_0'], step_frequency=fileLoad.step_frequency)

# plt.figure()
# input_pulse.plot_time(show=False)
#
# plt.figure()
# input_pulse.plot_frequency(show=False)

# plt.figure()
for i in range(len(fileLoad.get_matrix()[0])):
    vertical_pos = i * fileLoad.ds
    data = fileLoad.create_data_collector((vertical_pos, 3.9))
    data.fft()

    if not got_frequencies:
        frequency_axis = data.fft_frequencies
        got_frequencies = True
    # data.plot_frequency(show=False)
    data_collectors.append(data.fft_amplitude)


plt.figure()
plt.ylabel("Position (metres)")
plt.xlabel("Frequency (rad/s)")
im = plt.imshow(data_collectors, extent=[np.min(frequency_axis), np.max(frequency_axis), 4, 0], aspect='auto')
plt.show()
