import numpy as np
import matplotlib.pyplot as plt
import math

matrix = np.zeros((500, 500))

# add draw function here

# example drawing a square in the middle (from row 33 to row 67, column 33 to column 67)
def draw_rect(top_left, bottom_right, magnitude):
    row_start, col_start = top_left[0], top_left[1]
    row_end, col_end = bottom_right[0], bottom_right[1]

    matrix[row_start:row_end, col_start:col_end] += magnitude


def draw_lens(center, radius, thickness, magnitude):
    center_i, center_j = center[0], center[1]
    radius = radius    # for each curve on the convex lens (both curves are part of some circle)
    thickness = thickness  # of the entire lens
    displacer = radius - (thickness // 2)

    for j in range(-radius, -displacer):
        for i in range(-radius, 1 + radius):
            if i**2 + j**2 < radius**2:
                # print(center_i + i, center_j + j)
                matrix[center_i + i][center_j + j + displacer] += magnitude

    for j in range(displacer, radius + 1):
        for i in range(-radius, 1 + radius):
            if i**2 + j**2 < radius**2:
                matrix[center_i + i][center_j + j - displacer] += magnitude


def draw_quarter_circle(center, radius, rectangle_length, magnitude):
    center_i, center_j = center[0], center[1]
    for i in range(radius):
        for j in range(radius):
            if i**2 + j**2 < radius**2:
                matrix[center_i - i][center_j + j] += magnitude

    # left rectangle portion of waveguide
    draw_rect((center_i - radius, center_j - rectangle_length), center, magnitude)

    # down rectangle portion of waveguide
    draw_rect(center, (center_i + rectangle_length, center_j + radius), magnitude)


draw_lens((len(matrix)//4, len(matrix)//4), 200, 50, magnitude=1)
draw_rect((10, 10), (80, 50), magnitude=1)
draw_quarter_circle((len(matrix)//2, len(matrix)//2), 50, rectangle_length=100, magnitude=1)
draw_quarter_circle((len(matrix)//2, len(matrix)//2), 20, rectangle_length=100, magnitude=-1)


def main():
    fig, ax = plt.subplots()
    image = ax.imshow(matrix, vmax=1, vmin=0, aspect='auto')
    ax.set_aspect('equal', 'box')
    ax.set_title('Drawing')

    plt.show()
main()
