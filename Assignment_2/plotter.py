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
RESULTS_DIR = os.path.join(SCRIPT_DIR, "results")

# Add new problems here as needed, with their respective configurations
CONFIGS = [
    {
        "prefix": "prob1",
        "title": "Problem 1: Steady-State Temperature Distribution",
        "L": 0.3,
        "W": 0.4,
    },
    {
        "prefix": "prob2a",
        "title": "Problem 2 (Case A): Symmetry Boundary Conditions",
        "L": 0.3,
        "W": 0.4,
    },
    {
        "prefix": "prob2b",
        "title": "Problem 2 (Case B): Quarter Domain Symmetry",
        "L": 0.3,
        "W": 0.4,
        "mirror_quarter": True,
    },
]


def setup_environment():
    """Ensures the results directory exists."""
    if not os.path.exists(RESULTS_DIR):
        print(f"[*] Creating '{RESULTS_DIR}' directory...")
        os.makedirs(RESULTS_DIR)


def generate_plots(config):
    """Generates contour and line plots for a specific configuration."""
    csv_file = os.path.join(RESULTS_DIR, f"{config['prefix']}_results.csv")

    # If the Fortran code hasn't generated this specific file yet, skip it.
    if not os.path.exists(csv_file):
        return False
    print(f"[*] Loading data for {config['prefix']}...")
    df = pd.read_csv(csv_file)

    # QUARTER DOMAIN MIRRORING
    if config.get("mirror_quarter"):
        print("    -> Expanding quarter domain to full domain via symmetry...")

        grid = df.pivot(index="Y", columns="X", values="T")
        T_q = grid.values

        # Mirror Horizontally
        T_right = np.fliplr(T_q[:, :-1])
        T_bottom_half = np.concatenate((T_q, T_right), axis=1)

        # Mirror Vertically
        T_top_half = np.flipud(T_bottom_half[:-1, :])
        T_full = np.concatenate((T_bottom_half, T_top_half), axis=0)

        # Rebuild the X and Y coordinate arrays
        X_full = np.linspace(0.0, config["L"], T_full.shape[1])
        Y_full = np.linspace(0.0, config["W"], T_full.shape[0])
        X_mesh, Y_mesh = np.meshgrid(X_full, Y_full)

        # 5. Overwrite the original dataframe with our new massive dataset
        df = pd.DataFrame(
            {"X": X_mesh.flatten(), "Y": Y_mesh.flatten(), "T": T_full.flatten()}
        )

    # Generate contour plot
    print("[*] Generating 2D Temperature Contour...")

    grid = df.pivot(index="Y", columns="X", values="T")
    X = grid.columns.values
    Y = grid.index.values
    T = grid.values

    plt.figure(figsize=(8, 10))
    contour = plt.contourf(X, Y, T, levels=50, cmap="plasma")
    contour_lines = plt.contour(X, Y, T, levels=20, colors="black", linewidths=0.5)
    plt.clabel(contour_lines, inline=True, fontsize=8, fmt="%1.1f")
    plt.colorbar(contour, label="Temperature $(\\degree C)$")

    plt.title(config["title"], fontsize=14)
    plt.xlabel("X Coordinate (m)", fontsize=12)
    plt.ylabel("Y Coordinate (m)", fontsize=12)
    plt.axis("equal")

    save_path = os.path.join(RESULTS_DIR, f"{config['prefix']}_contour.png")
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"    -> Saved to {save_path}")
    plt.close()

    # Generate line plots
    print("[*] Generating 0.05m Interval Line Plots...")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    max_dim = max(config["L"], config["W"])
    target_intervals = np.arange(0.0, max_dim + 0.05, 0.05)

    # Plot T vs Y at specific X intervals
    for x_val in target_intervals:
        if x_val <= config["L"]:
            slice_data = df[np.isclose(df["X"], x_val, atol=1e-5)]
            if not slice_data.empty:
                ax1.plot(slice_data["Y"], slice_data["T"], label=f"X = {x_val:.2f} m")

    ax1.set_title("Temperature Profile along Y (Constant X slices)")
    ax1.set_xlabel("Y Coordinate (m)")
    ax1.set_ylabel("Temperature $(\\degree C)$")
    ax1.grid(True, linestyle="--", alpha=0.7)
    ax1.legend(loc="upper right", fontsize=8)

    # Plot T vs X at specific Y intervals
    for y_val in target_intervals:
        if y_val <= config["W"]:
            slice_data = df[np.isclose(df["Y"], y_val, atol=1e-5)]
            if not slice_data.empty:
                ax2.plot(slice_data["X"], slice_data["T"], label=f"Y = {y_val:.2f} m")

    ax2.set_title("Temperature Profile along X (Constant Y slices)")
    ax2.set_xlabel("X Coordinate (m)")
    ax2.set_ylabel("Temperature $(\\degree C)$")
    ax2.grid(True, linestyle="--", alpha=0.7)
    ax2.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=4, fontsize=8)

    line_path = os.path.join(RESULTS_DIR, f"{config['prefix']}_line_plots.png")
    plt.savefig(line_path, dpi=300, bbox_inches="tight")
    print(f"    -> Saved to {line_path}")
    plt.close()

    return True


if __name__ == "__main__":
    print("=========================================================")
    print("                  CFD Post-Processing                    ")
    print("=========================================================")

    # Check directories and get data path
    setup_environment()

    # Keep track of how many files we successfully processed to determine if we need to show an error message at the end
    files_processed = 0

    for config in CONFIGS:
        if generate_plots(config):
            files_processed += 1

    if files_processed == 0:
        print("\n[!] No data files found in the results directory.")
        print("---------------------------------------------------------")
        print(" ACTION REQUIRED:")
        print(" 1. Run your compiled Fortran executable.")
        print(" 2. Run this Python script again to generate the plots.")
        print("---------------------------------------------------------\n")
        sys.exit(0)

    print("=========================================================")
    print("[SUCCESS] All post-processing complete.")
