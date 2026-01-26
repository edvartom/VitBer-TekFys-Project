import numpy as np
import matplotlib.pyplot as plt
from funcs import f_1c, RK32

x_init = 0
x_end = 2*np.pi
y_init = np.array([0,2],dtype=np.float64)
h0 = 0.005
tol = 1e-7
alpha = 0.8

x_arr, y_arr, h_arr, N = RK32(x_init, x_end, y_init, f_1c, h0, tol, alpha)


fig_1c, ax_1c = plt.subplots(2,1)
ax_1c[0].plot(x_arr, np.sin(2*x_arr), color="cyan",label='y(x) analytic')
ax_1c[0].plot(x_arr, 2*np.cos(2*x_arr), color="cyan",label="y'(x) analytic")
ax_1c[0].plot(x_arr,y_arr[:,0],color="tab:red",linestyle="--",label='y')
ax_1c[0].plot(x_arr,y_arr[:,1],color="orange",linestyle="--",label="y'(x)")
ax_1c[0].set_title("Solution of differential equation\ny''(x)=-4sin(2x)")
ax_1c[0].set_xlabel('x')
ax_1c[0].set_ylabel('y(x)')
ax_1c[0].legend()
ax_1c[1].plot(x_arr,h_arr)
ax_1c[1].set_title("Steplengths used")
ax_1c[1].set_xlabel('x')
ax_1c[1].set_ylabel('Steplength, h')

plt.tight_layout()
plt.show()