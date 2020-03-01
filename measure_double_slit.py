from analyser import FileLoader, DataCollector
import matplotlib.pyplot as plt
import numpy as np

# fileLoad = FileLoader("large_double_slit")
fileLoad = FileLoader("far_spaced_double_slit")

data_collectors = []

got_frequencies = False
frequency_x = None

for i in range(len(fileLoad.get_matrix()[0])):
    data = DataCollector(fileLoad.get_matrix(), i, -5)
    data.collect_all()
    data.fft()

    if not got_frequencies:
        frequency_x = data.fft_frequencies
        got_frequencies = True
    #
    # data.plot_frequency()
    # data.plot_time()

    data_collectors.append(data.fft_amplitude)


print(len(data_collectors))
plt.figure()
plt.title("2 slit diffraction")

plt.ylabel("Position (metres)")
plt.xlabel("Frequency (per timestep)")
# for i in range(len(data_collectors)):
#     plt.plot(data_collectors[i])
im = plt.imshow(data_collectors, extent=[frequency_x[0], frequency_x[-1], 4, 0], aspect='auto')
plt.show()

