import numpy as np
import matplotlib.pyplot as plt

# Parameters for constantan block
k = 23  # Thermal conductivity for constantan in W/m-K
h = 1000  # Convection heat transfer coefficient in W/m^2-K, an assumed value
q_top = 8000  # Heat flux from the heater in W
block_length = 0.5  # Block length in m
block_height = 0.3  # Block height in m
T_bottom = T_sides = 0  # Temperature at the bottom and sides in °C

# Grid parameters
nx, ny = 7, 5  # Number of blocks in x and y directions
dx = block_length / (nx - 1)
dy = block_height / (ny - 1)

# Initialize temperature grid
T = np.full((ny, nx), T_bottom)

# Apply boundary conditions
T[-1, :] = q_top / (k * dy)  # Top boundary condition from the heat flux
T[:, 0] = T[:, -1] = T_bottom  # Side boundary conditions

# Iterative parameters
max_iter = 10000
tolerance = 1e-6
error = np.inf
iterations = 0

# Iterative solver for steady-state
while error > tolerance and iterations < max_iter:
    iterations += 1
    T_old = T.copy()

    # Update temperature for each internal node
    for i in range(1, ny-1):
        for j in range(1, nx-1):
            T[i, j] = (T[i+1, j] + T[i-1, j] + T[i, j+1] + T[i, j-1]) / 4.0

    # Compute the maximum change in temperature
    error = np.abs(T - T_old).max()

# Output results
print(f'Steady-state solution reached after {iterations} iterations with a maximum change of {error:.6f} degrees.')
print(f'Maximum temperature in the grid is {T.max()} °C')

# Plot the results
plt.imshow(T, cmap='jet', interpolation='nearest', origin='lower', extent=[0, block_length, 0, block_height], vmin=np.min(T), vmax=np.max(T))
plt.colorbar(label='Temperature (°C)')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.title('Steady-state Temperature Distribution on a 7x5 Grid')
plt.show()
