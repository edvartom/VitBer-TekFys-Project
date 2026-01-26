from funcs import *
import numpy as np
import matplotlib.pyplot as plt
import time

# Defining params
h0 = 0.005
tol = 1e-7
alpha = 0.8
x_init = 0              # Starting position
x_end = 12              # Ending position

b0 = 0                  # Initial guess 0
b1 = 8                  # Initial guess 1

start_g = time.perf_counter() # Measuring time for later use

x_arr_g, y_arr_g, h_arr_g, N_g = BVP_solver(secant_method, RK32, b0, b1, x_init, x_end, f_1g, h0, tol, alpha)

stop_g = time.perf_counter()

runtime_g = start_g - stop_g

plt.plot(x_arr_g, y_arr_g[:,0], label="$y(x)$")
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.title("Numerical solution of BVP")
plt.legend()
plt.show()