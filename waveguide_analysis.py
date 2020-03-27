from analyser import FileLoader
from solver import Pulse
from solver import MaterialArray
import matplotlib.pyplot as plt

fileLoad = FileLoader("waveguide")
material = MaterialArray(fileLoad.size, fileLoad.length_x, fileLoad.length_y, eps_array=fileLoad.material)
fileLoad.play(1)
# material.plot()

data_output = fileLoad.create_data_collector((2.85, 2.6))
data_input = fileLoad.create_data_collector(fileLoad.constants['pulses'][0]['location'])

plt.figure()
data_input.plot_time(show=False)
data_output.plot_time(show=False)

plt.figure()
data_output.plot_frequency(show=False)
data_input.plot_frequency(show=False)

plt.figure()
transfer_func = data_input.fft_amplitude * data_output.fft_amplitude

plt.plot(data_output.fft_frequencies, transfer_func)
plt.show()