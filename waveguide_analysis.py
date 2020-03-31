from analyser import FileLoader
from solver import Pulse
from solver import MaterialArray
import matplotlib.pyplot as plt

fileLoad = FileLoader("slab_waveguide")
material = MaterialArray(fileLoad.size, fileLoad.length_x, fileLoad.length_y, eps_array=fileLoad.material)
material.plot()
fileLoad.play(1)

data_output = fileLoad.create_data_collector((3.5, 3.5))
data_input = fileLoad.create_data_collector((fileLoad.constants['pulses'][0]['location'][0], 0.1 + fileLoad.constants['pulses'][0]['location'][1]))

plt.figure()
data_input.plot_time(show=False)
data_output.plot_time(show=False)

plt.figure()
data_output.plot_frequency(show=False)
data_input.plot_frequency(show=False)
# data_output.fft()
# data_input.fft()
#
# plt.figure()
# transfer_mask = data_output.fft_frequencies < 5 * 10**9
# transfer_func = data_output.fft_amplitude / data_input.fft_amplitude
# transfer_func[transfer_mask] = 0
# plt.plot(data_output.fft_frequencies, transfer_func)
#
# plt.figure()
# multiply_func = data_output.fft_amplitude * data_input.fft_amplitude
# plt.plot(data_output.fft_frequencies, multiply_func)


plt.show()
