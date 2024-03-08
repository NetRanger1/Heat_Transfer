#This is the code for 5-90
import numpy as np
import matplotlib.pyplot as plt

print("2D heat equation solver")

# Define the physical dimensions of the plate (in meters)
x_length = 0.09  # 9 cm
y_length = 0.045  # 4.5 cm

# Thermal conductivity values for different materials (W/m-K)
kappa_values = {'Cu': 401, 'Ni': 90.7, 'Fe': 15.1}

# Heat transfer coefficient (W/m^2-K)
h = 125

# Spatial step size (in meters)
delta_x = 0.015  # 1.5 cm

# Calculate the number of nodes in each direction, adding 1 to include both boundaries
x_norm = round(x_length / delta_x) + 1
y_norm = round(y_length / delta_x) + 1
print("nodes in x: ", x_norm)
print("nodes in y: ", y_norm)

# Initialize the temperature grid with an initial condition and boundary conditions
T = np.empty((y_norm, x_norm))  # Temperature array
T_initial = 293  # Initial temperature inside the grid (in Kelvin)
T_amb = 288  # Ambient temperature outside the grid (in Kelvin)

# Set the initial condition across the grid
T.fill(T_initial)

# Set fixed temperatures along the bottom boundary of the grid (in Kelvin)
T[:1, :] = np.array([600, 700, 800, 900, 800, 700, 600]) + 273

# Iteration control variables
max_T_change = 1  # Maximum temperature change in an iteration step
iter_to_solve = 0  # Counter for the number of iterations

# Main loop to iteratively solve the heat equation until convergence
while max_T_change >= 0.001:
    iter_to_solve += 1
    temp_T = np.copy(T)  # Temporary copy for comparison

    # Iterate over the interior nodes to update temperatures
    for i in range(1, y_norm - 1):
        for j in range(1, x_norm - 1):
            # Determine material based on position
            material = 'Cu' if 1 < j <= 3 else 'Ni' if 3 < j <= 5 else 'Fe'
            kappa = kappa_values[material]  # Thermal conductivity

            # Top boundary condition using convective heat transfer formula
            if i == y_norm - 1:
                T[i][j] = ((h * delta_x / kappa * T_amb) + T[i - 1][j] + T[i][j + 1] + T[i][j - 1]) / (
                            h * delta_x / kappa + 2)
            else:
                # Update temperature using the discrete Laplace equation (average of neighbors)
                T[i][j] = (T[i + 1][j] + T[i - 1][j] + T[i][j + 1] + T[i][j - 1]) / 4

    # Calculate the maximum temperature change in this iteration
    max_T_change = np.max(np.abs(T - temp_T))
    print("Max delta T: ", max_T_change)

print("Maximum Temperature : ", np.max(T))
print(f"Solution took {iter_to_solve:.0f} iterations to complete")

# Adjust temperatures to Celsius for plotting
T -= 273

# Plotting the temperature distribution
plt.clf()
plt.title("Temperature (C)")
plt.xlabel("x")
plt.ylabel("y")
x_range = np.arange(0, x_length + delta_x, delta_x)
y_range = np.arange(0, y_length + delta_x, delta_x)
plt.pcolormesh(x_range, y_range, T, cmap=plt.cm.jet, vmin=np.min(T), vmax=np.max(T))
plt.colorbar()
plt.savefig('5-90_SS.png')
plt.show()
