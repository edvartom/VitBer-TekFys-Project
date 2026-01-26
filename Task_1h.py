from funcs import *
from Task_1g import x_arr_g, y_arr_g, runtime_g
import numpy as np
import time
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt

# Defining constants and variables
x_init = 0
x_end = 12
x = np.linspace(x_init, x_end, 10)      # Initial x-mesh
y = np.zeros((2, x.size))               # Boundaries for y, 2D array, x rows.
start_h = time.perf_counter()           # Timing scipy BVP solver

# Solving BVP with scipy
def bc(ya, yb): # Ensures y(0)=0, y(12)=0
    return np.array([ya[0], yb[0]])

y_solution = solve_bvp(f_1g, bc, x, y) # Scipy solver
y_h = y_solution.sol(x_arr_g)[0] # Fetches solution for y(x). Allows integrated interpolation for smoother plots.

stop_h = time.perf_counter() 

diff_arr = np.abs(y_h - y_arr_g[:,0]) # Absolute differences between 1g and 1h solutions

# Plotting:
plt.plot(x_arr_g, y_h, 'g', label="Scipy")
plt.plot(x_arr_g, y_arr_g[:,0], 'r--', label="RK32")
plt.title("y(x) solved by Scipy and RK3")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.grid()
plt.show()

plt.plot(x_arr_g, diff_arr, label="Difference as function of x")
plt.title("Error difference between scipy BVP solver, and our BVP solver")
plt.xlabel("x")
plt.ylabel("Difference in y")
plt.legend()
plt.show()

print(f"Scipy's BVP was faster than our solver by {start_h - stop_h - runtime_g} seconds")