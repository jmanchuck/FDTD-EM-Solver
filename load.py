import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

h_list = np.load('h_list.npy')

def update(i):
    matrice.set_array(h_list[i])

fig, ax = plt.subplots()
matrice = ax.imshow(h_list[0], vmax=1 * 10**8, vmin= - 1 * 10**8)
plt.colorbar(matrice)

ani = animation.FuncAnimation(fig, update, frames=len(h_list), interval=100, repeat=True)
plt.show()
