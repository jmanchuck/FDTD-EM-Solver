from analyser import FileLoader, DataCollector
import matplotlib.pyplot as plt
import numpy as np

fileLoad = FileLoader("double_slit")
# fileLoad.play(1)

data_collectors = []

get_freqs = False
frequency_x = None

for i in range(len(fileLoad.get_matrix()[0])):
    data = DataCollector(fileLoad.get_matrix(), i, -10)
    data.collect_all()
    data.fft()

    if not get_freqs:
        frequency_x = data.fft_frequencies
        get_freqs = True

    data_collectors.append(data.fft_amplitude)


# print(len(data_collectors))
plt.figure()
for i in range(len(data_collectors) // 10):
    plt.plot(data_collectors[10 * i])
# im = plt.imshow(data_collectors)
plt.show()

