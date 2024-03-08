#This is the code for 5-78
import numpy as np
import matplotlib.pyplot as plt

# Setting up parameters for the heat conduction problem in a constantan block
k = 23  # Thermal conductivity of constantan in W/m-K
h = 1000  # Convection heat transfer coefficient in W/m^2-K (assumed value)
q_top = 8000  # Heat flux from the heater at the top surface in W
block_length = 0.5  # Length of the block in meters
block_height = 0.3  # Height of the block in meters
T_bottom = T_sides = 0  # Initial temperature at the bottom and sides in °C, assuming a cold start

# Defining the discretization grid
nx, ny = 7, 5  # Number of divisions in x (length) and y (height) directions
dx = block_length / (nx - 1)  # Distance between nodes in x-direction
dy = block_height / (ny - 1)  # Distance between nodes in y-direction

# Initializing the temperature grid with the bottom and sides temperature
T = np.full((ny, nx), T_bottom)  # Initialize the temperature array with the boundary condition value

# Applying boundary conditions
T[-1, :] = q_top / (k * dy)  # Top boundary condition derived from the heat flux
T[:, 0] = T[:, -1] = T_bottom  # Side boundary conditions (temperatures)

# Iteration setup for solving the system
max_iter = 1000  # Maximum number of iterations
tolerance = 1e-6  # Convergence criterion
error = np.inf  # Initialize error to a large value
iterations = 0  # Iteration counter

# Iteratively solving for the temperature distribution until steady-state is reached
while error > tolerance and iterations < max_iter:
    iterations += 1  # Increment iteration counter
    T_old = T.copy()  # Make a copy of the temperature array for comparison

    # Updating the temperature for internal nodes using the finite difference method
    for i in range(1, ny-1):
        for j in range(1, nx-1):
            T[i, j] = (T[i+1, j] + T[i-1, j] + T[i, j+1] + T[i, j-1]) / 4.0

    # Compute the maximum change in temperature across the grid to check for convergence
    error = np.abs(T - T_old).max()

# Output the results of the iterative process
print(f'Steady-state solution reached after {iterations} iterations with a maximum change of {error:.6f} degrees.')
print(f'Maximum temperature in the grid is {T.max()} °C')

# Plotting the steady-state temperature distribution
plt.imshow(T, cmap='jet', interpolation='nearest', origin='lower',
           extent=[0, block_length, 0, block_height],
           vmin=np.min(T), vmax=np.max(T))  # Setting the plot's color scale based on the temperature range
plt.colorbar(label='Temperature (°C)')  # Color bar to indicate temperature values
plt.xlabel('x (m)')  # Label for the x-axis
plt.ylabel('y (m)')  # Label for the y-axis
plt.title('Steady-state Temperature Distribution on a 7x5 Grid')  # Title for the plot
plt.savefig('5-78_SS.png')  # Save the plot as a PNG file
plt.show()  # Display the plot
