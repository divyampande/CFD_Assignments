import numpy as np
import matplotlib.pyplot as plt
import os

# Input parameters
l = 1.0  # Length of the domain
alpha = 1.0  # Diffusivity
nodes = 11  # Number of grid points
dx = l / (nodes - 1)  # Grid spacing
dt = 0.1  # Time step size
gamma = alpha * dt / dx**2  # Stability parameter
total_time = 10.0  # Total simulation time

# Check for numerical stability
if gamma > 0.5:
    print(
        f"WARNING: Stability condition (gamma <= 0.5) violated. Current gamma = {gamma:.3f}"
    )
    print("Consider reducing 'dt' or increasing 'dx' to ensure numerical stability.")
num_steps = int(total_time / dt)  # Number of time steps

time_plot = np.array(
    [0, 0.002, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0, 2.0, 5.0, 10.0]
)  # Array to store the time steps which will be plotted

# Spatial coordinates
x = np.linspace(0, l, nodes)


# Initial and boundary conditions
T = np.zeros(nodes)  # Initial temperature distribution
T[0] = 1.0  # Boundary condition at the left end
T[-1] = 0.0  # Boundary condition at the right end

# Array to store temperature distribution
T_history = np.zeros((num_steps + 1, nodes))
T_history[0] = T.copy()  # Store initial condition
T_history[:, 0] = T[0].copy()  # Store left boundary condition
T_history[:, -1] = T[-1].copy()  # Store right boundary condition

# Time-stepping loop
for n in range(1, num_steps + 1):
    T_old = T.copy()
    T[1:-1] = T_old[1:-1] + gamma * (T_old[2:] - 2 * T_old[1:-1] + T_old[:-2])
    T_history[n] = T.copy()

# Since the left boundary condition is applied for t>0 but not t=0, we need to update the first row of T_history.
# The left boundary is still 0 at t=0.
T_history[0, 0] = 0.0


def plot_temperature_distribution(
    x_coords, T_history, time_points, dt, gamma, output_filename
):
    plt.style.use("seaborn-v0_8-whitegrid")
    fig, ax = plt.subplots(figsize=(12, 7))

    # Use a colormap for visually distinct line colors
    colors = plt.cm.viridis(np.linspace(0, 1, len(time_points)))

    for i, t in enumerate(time_points):
        # find the closest time step index to avoid floating point issues
        step_index = int(round(t / dt))

        if step_index < T_history.shape[0]:
            ax.plot(
                x_coords,
                T_history[step_index],
                marker="o",
                markersize=4,
                linestyle="-",
                color=colors[i],
                label=f"t = {t:.3f} s",
            )

    ax.set_xlabel("Position (m)", fontsize=12)
    ax.set_ylabel("Temperature", fontsize=12)
    ax.set_title(
        f"1D Heat Conduction (FTCS) | $\Delta t$ = {dt} s | $\gamma$ = {gamma:.3f}",
        fontsize=14,
        fontweight="bold",
    )
    ax.legend(title="Time", fontsize=10)
    plt.savefig(output_filename, dpi=300, bbox_inches="tight")
    print(f"Plot saved successfully as: {output_filename}")
    plt.show()


figname = f"temperature_distribution_at_dt_{dt}.png"
script_dir = os.path.dirname(__file__)
plot_path = os.path.join(script_dir, figname)
plot_temperature_distribution(x, T_history, time_plot, dt, gamma, plot_path)
