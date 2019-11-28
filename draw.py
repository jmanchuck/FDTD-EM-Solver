import numpy as np
import matplotlib.pyplot as plt
import math

matrix = np.zeros((500, 500))

################################################
# add draw function here

# example drawing a square in the middle (from row 33 to row 67, column 33 to column 67)
# matrix[33:67, 33:67] = 1

center_i, center_j = (len(matrix) // 2, len(matrix) // 2)
radius = 100    # for each curve on the convex lens (both curves are part of some circle)
thickness = 50  # of the entire lens
displacer = radius - (thickness // 2)

for j in range(-radius, -displacer):
    for i in range(-radius, 1 + radius):
        if i**2 + j**2 < radius**2:
            # print(center_i + i, center_j + j)
            matrix[center_i + i][center_j + j + displacer] = 1

for j in range(displacer, radius + 1):
    for i in range(-radius, 1 + radius):
        if i**2 + j**2 < radius**2:
            matrix[center_i + i][center_j + j - displacer] = 1

################################################

def main():
    fig, ax = plt.subplots()
    image = ax.imshow(matrix, vmax=1, vmin=0, aspect='auto')
    ax.set_aspect('equal', 'box')
    ax.set_title('Drawing')

    plt.show()
main()
