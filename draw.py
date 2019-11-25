import numpy as np
import matplotlib.pyplot as plt

matrix = np.zeros((100, 100))

################################################
# add draw function here

# example drawing a square in the middle (from row 33 to row 67, column 33 to column 67
matrix[33:67, 33:67] = 1


################################################

def main():
    fig, ax = plt.subplots()
    image = ax.imshow(matrix, vmax=1, vmin=0, aspect='auto')
    ax.set_aspect('equal', 'box')
    ax.set_title('Drawing')

    plt.show()
main()
