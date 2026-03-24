# ======================================================================
# File: plotter.py
# Author: Divyam Pandey (AE25M021)
# Description: Reads CFD output data from Fortran and generates
#              2D temperature contours and specific line profiles.
# ======================================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_temperature_contour(csv_filename, title, L, W):
    # TODO: Load CSV using pandas/numpy
    # TODO: Create X, Y meshgrid based on L and W
    # TODO: Use plt.contourf() to plot the heat map
    pass


def plot_line_profiles(csv_filename, L, W):
    # TODO: Extract specific rows/columns corresponding to 0.05m steps
    # TODO: Generate line plots as requested in Problem 1 [cite: 30]
    pass


if __name__ == "__main__":
    print("Generating CFD Plots...")

    plt.show()
    print("Done!")
