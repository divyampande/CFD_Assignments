# ======================================================================
# File: plotter.py
# Author: Divyam Pandey (AE25M021)
# Description: Reads CFD output data from Fortran and generates
#              2D temperature contours and specific line profiles.
# ======================================================================
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Get the absolute directory path where this specific python file lives
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Lock the results directory to that exact location
RESULTS_DIR = os.path.join(SCRIPT_DIR, "results")


def setup_environment():
    """Checks for the results directory and required data files."""
    if not os.path.exists(RESULTS_DIR):
        print(f"[*] Creating '{RESULTS_DIR}' directory...")
        os.makedirs(RESULTS_DIR)

    prob1_csv = os.path.join(RESULTS_DIR, "prob1_results.csv")

    if not os.path.exists(prob1_csv):
        print(f"\n[!] Data file '{prob1_csv}' not found.")
        print("---------------------------------------------------------")
        print(" ACTION REQUIRED:")
        print(" 1. The 'results' directory is ready.")
        print(" 2. Please run your compiled Fortran executable now.")
        print(" 3. Run this Python script again to generate the plots.")
        print("---------------------------------------------------------\n")
        sys.exit(0)
    return prob1_csv


def plot_temperature_contour(df):
    """Generates a 2D color contour map of the steady-state temperature."""
    print("[*] Generating 2D Temperature Contour...")

    grid = df.pivot(index="Y", columns="X", values="T")
    X = grid.columns.values
    Y = grid.index.values
    T = grid.values

    plt.figure(figsize=(8, 10))
    contour = plt.contourf(X, Y, T, levels=100, cmap="inferno")
    plt.colorbar(contour, label="Temperature (°C)")

    plt.title("Problem 1: Steady-State Temperature Distribution", fontsize=14)
    plt.xlabel("X Coordinate (m)", fontsize=12)
    plt.ylabel("Y Coordinate (m)", fontsize=12)
    plt.axis("equal")

    save_path = os.path.join(RESULTS_DIR, "prob1_contour.png")
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"    -> Saved to {save_path}")


def plot_line_profiles(df):
    """Plots temperature at 0.05m intervals along X and Y axes."""
    print("[*] Generating 0.05m Interval Line Plots...")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    target_intervals = np.arange(0.0, 0.45, 0.05)

    # Plot T vs Y at specific X intervals
    for x_val in target_intervals:
        if x_val <= 0.3:  # L is max 0.3m
            slice_data = df[np.isclose(df["X"], x_val, atol=1e-5)]
            if not slice_data.empty:
                ax1.plot(slice_data["Y"], slice_data["T"], label=f"X = {x_val:.2f} m")

    ax1.set_title("Temperature Profile along Y (Constant X slices)")
    ax1.set_xlabel("Y Coordinate (m)")
    ax1.set_ylabel("Temperature (°C)")
    ax1.grid(True, linestyle="--", alpha=0.7)
    ax1.legend(loc="upper right", fontsize=8)

    # Plot T vs X at specific Y intervals
    for y_val in target_intervals:
        if y_val <= 0.4:  # W is max 0.4m
            slice_data = df[np.isclose(df["Y"], y_val, atol=1e-5)]
            if not slice_data.empty:
                ax2.plot(slice_data["X"], slice_data["T"], label=f"Y = {y_val:.2f} m")

    ax2.set_title("Temperature Profile along X (Constant Y slices)")
    ax2.set_xlabel("X Coordinate (m)")
    ax2.set_ylabel("Temperature (°C)")
    ax2.grid(True, linestyle="--", alpha=0.7)
    ax2.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=4, fontsize=8)

    save_path = os.path.join(RESULTS_DIR, "prob1_line_plots.png")
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"    -> Saved to {save_path}")


if __name__ == "__main__":
    print("=========================================================")
    print("                  CFD Post-Processing                    ")
    print("=========================================================")

    # Check directories and get data path
    data_file = setup_environment()

    print(f"[*] Loading data from {data_file}...")
    df = pd.read_csv(data_file)

    plot_temperature_contour(df)
    plot_line_profiles(df)

    print("=========================================================")
    print("[SUCCESS] All post-processing complete.")

    # Display the plots
    # plt.show()
