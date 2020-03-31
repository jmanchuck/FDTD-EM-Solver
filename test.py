from analyser import FileLoader as FL

fL = FL("test")

# fL.play()

inpt = fL.create_data_collector((0.7, 0.3))
outpt = fL.create_data_collector((3.5, 3))

inpt.fft(), outpt.fft()

inpt.plot_time(show=False)
outpt.plot_time(show=True)