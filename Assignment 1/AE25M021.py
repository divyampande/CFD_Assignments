import numpy as np
import matplotlib.pyplot as plt

# Input parameters
l = 1.0  # Length of the domain
alpha = 1.0  # Diffusivity
nodes = 11  # Number of grid points
dx = l / (nodes - 1)  # Grid spacing
dt = 0.001  # Time step
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
    [0, 0.1, 0.5, 1, 2, 5, 10]
)  # Array to store the time steps which will be plotted


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
    # Vectorized update for inner nodes
    T[1:-1] = T_old[1:-1] + gamma * (T_old[2:] - 2 * T_old[1:-1] + T_old[:-2])
    T_history[n] = T.copy()

# Final plot of temperature distribution at different time steps
plt.figure(figsize=(10, 6))
for t in time_plot:
    plt.plot(
        np.linspace(0, l, nodes), T_history[int(t / dt)], label=f"Time = {t:.3f} s"
    )
plt.legend()
plt.xlabel("Position (m)")
plt.ylabel("Temperature (°C)")
plt.title("Temperature Distribution at Different Time Steps")
plt.grid()
plt.show()
