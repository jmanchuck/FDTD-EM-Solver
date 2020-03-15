from analyser import FileLoader, DataCollector
import matplotlib.pyplot as plt
import numpy as np

# fileLoad = FileLoader("large_double_slit")
fileLoad = FileLoader("far_spaced_double_slit")

data_collectors = []

got_frequencies = False
frequency_x = None

for i in range(len(fileLoad.get_matrix()[0])):
    data = DataCollector(fileLoad.get_matrix(), i, -5, fileLoad.constants['dt'], fileLoad.constants['end_time'], fileLoad.constants['step_frequency'])
    data.collect_all()
    data.fft()

    if not got_frequencies:
        frequency_x = data.fft_frequencies
        got_frequencies = True

    data_collectors.append(data.fft_amplitude)
plt.figure()
plt.title("2 slit diffraction")

plt.ylabel("Position (metres)")
plt.xlabel("Frequency (per timestep)")
im = plt.imshow(data_collectors, extent=[frequency_x[0], fileLoad.constants["pulses"][0]["sigma_w"], 9, 0], aspect='auto')

plt.show()
