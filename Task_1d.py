from funcs import *
from scipy.interpolate import interp1d

# Parameters for errorplot and timestep-plot:
x_init = 0
x_end = 2*np.pi
y_init = np.array([0,2],dtype=np.float64)
h0 = 0.0003
tol_arr = np.linspace(1e-7,1e-1,10)
tol = 1e-7
alpha = 0.8
alpha_arr = np.linspace(0.1,1.0,10)

# Values errorplot:
error_arr = error_afo_tol(x_init, x_end, y_init, f_1c,h0, 
                          tol_arr, alpha, y_anl_1)
# Values timestep-plot:
N_timesteps_arr = N_timesteps_afo_alpha(x_init, x_end, y_init, f_1c,h0, 
                                        tol, alpha_arr)

# Plotting errorplot
fig_1d,ax_1d=plt.subplots(2,1)
ax_1d[0].scatter(tol_arr, error_arr)
ax_1d[0].set_title('Error as function of time tol (error tolerance)')
ax_1d[0].set_xlabel('tolerance, tol')
ax_1d[0].set_ylabel('maximum value of error')
ax_1d[0].set_xticks(np.linspace(tol_arr[0], tol_arr[-1], 7))
# Plotting timestep-plot
ax_1d[1].scatter(alpha_arr, N_timesteps_arr)
ax_1d[1].set_title('Number of timesteps (both used and discarded) \na.f.o. the pessimist factor, alpha')
ax_1d[1].set_xlabel('alpha')
ax_1d[1].set_ylabel('timesteps, N')
ax_1d[1].set_xticks(alpha_arr)
#
plt.tight_layout()
plt.show()


"""

# Simple Visualisation to determine linearity
m = (error_arr[-1] - error_arr[0]) / (tol_arr[-1] - tol_arr[0])
x = np.linspace(tol_arr[0], tol_arr[-1], 20)
y = m * x + error_arr[0] 

plt.plot(x, y)
plt.plot(tol_arr, error_arr, 'bo')
plt.show()
"""