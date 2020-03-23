from analyser import FileLoader
import matplotlib.pyplot as plt
import numpy as np

# fileLoad = FileLoader("far_spaced_double_slit")
fileLoad = FileLoader("d_slit_test")

data_collectors = []

got_frequencies = False
frequency_axis = None

for i in range(len(fileLoad.get_matrix()[0])):
    vertical_pos = i * fileLoad.ds
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
