from analyser import FileLoader
import matplotlib.pyplot as plt
import numpy as np

fileLoad = FileLoader("d_slit_test")

data = fileLoad.create_data_collector((2, 0.2))
data.collect_all()
data.fft()
data.plot_frequency()
print(len(data.fft_frequencies))
plt.show()