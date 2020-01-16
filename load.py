import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

h_list = np.load('double_slit.npy')

max = np.max(h_list)

def update(i):
    matrice.set_array(h_list[i])

fig, ax = plt.subplots()
matrice = ax.imshow(h_list[0], vmax=0.7*max, vmin=-0.7*max)
plt.colorbar(matrice)

ani = animation.FuncAnimation(fig, update, frames=len(h_list), interval=150, repeat=True)
plt.show()
