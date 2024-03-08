import numpy as np
import matplotlib.pyplot as plt

print("2D heat equation solver")
x_length = 0.09  # 9 cm by 4.5 cm
y_length = 0.045
kappa_values = {'Cu': 401, 'Ni': 90.7, 'Fe': 15.1}  # W/m-K
h = 125  # W/m^2-K
delta_x = 0.015  # 1.5 cm
x_norm = round(x_length / delta_x) + 1  # node on each end need to add 1
y_norm = round(y_length / delta_x) + 1
print("nodes in x: ", x_norm)
print("nodes in y: ", y_norm)

# Initialize solution: the grid of u(k, i, j)
T = np.empty((y_norm, x_norm))
T_initial = 293  # Initial condition everywhere inside the grid
T_amb = 288  # Temperature in K

# Set the initial condition
T.fill(T_initial)
T[:1, :] = np.array([600, 700, 800, 900, 800, 700, 600]) + 273  # fixed T on bottom

max_T_change = 1
iter_to_solve = 0

while max_T_change >= 0.001:
    iter_to_solve += 1
    temp_T = np.copy(T)
    for i in range(1, y_norm - 1):
        for j in range(1, x_norm - 1):
            material = 'Cu' if 1 < j <= 3 else 'Ni' if 3 < j <= 5 else 'Fe'
            kappa = kappa_values[material]
            if i == y_norm - 1:  # top boundary condition
                T[i][j] = ((h * delta_x / kappa * T_amb) + T[i - 1][j] + T[i][j + 1] + T[i][j - 1]) / (
                        h * delta_x / kappa + 2)
            else:  # bulk
                T[i][j] = (T[i + 1][j] + T[i - 1][j] + T[i][j + 1] + T[i][j - 1]) / 4

    max_T_change = np.max(np.abs(T - temp_T))
    print("Max delta T: ", max_T_change)

print("Maximum Temperature : ", np.max(T))
print(f"Solution took {iter_to_solve:.0f} iterations to complete")

# convert to K from C
T -= 273

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
