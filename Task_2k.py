from funcs import *


delta=0.01
m=101
l_arr=np.array([0.5,1,2]) # l values to test
epsilon_arr=np.linspace(0,2,m)[::-1] # Epsilon array going from 2 to 0
y_init=np.zeros((32,m)) # Initial guess (for eps=2)

# Plotting
fig_2k, ax_2k=plt.subplots(nrows=3, ncols=1, figsize=(12,6))
fig_2k.suptitle("Density of states for different energies at $x=l/2$ for $l=\{0.5,1,2\}$")
x_init=np.linspace(0,l_arr[0],m) # Position
DOS_at_mid05,sol_arr05=DOS_diff_eps(epsilon_arr, x_init, y_init,l_arr[0],m)

ax_2k[0].set_title("$l=0.5$")
ax_2k[0].plot(epsilon_arr,DOS_at_mid05,color="red")
ax_2k[0].set_ylabel("$DOS, D/D_0$")
ax_2k[0].set_xlabel("epsilon")
ax_2k[0].grid()

x_init=np.linspace(0,l_arr[1],m) # Position
DOS_at_mid1,sol_arr1=DOS_diff_eps(epsilon_arr, x_init, y_init,l_arr[1],m)
ax_2k[1].plot(epsilon_arr,DOS_at_mid1,color="red")
ax_2k[1].set_ylabel("$DOS, D/D_0$")
ax_2k[1].set_xlabel("epsilon")
ax_2k[1].grid()

x_init=np.linspace(0,l_arr[2],m) # Position
DOS_at_mid2,sol_arr2=DOS_diff_eps(epsilon_arr, x_init, y_init,l_arr[2],m)
ax_2k[2].plot(epsilon_arr,DOS_at_mid2,color="red")
ax_2k[2].set_ylabel("$DOS, D/D_0$")
ax_2k[2].set_xlabel("epsilon")
ax_2k[2].grid()

ax_2k[0].tight_layout()
ax_2k[1].tight_layout()
ax_2k[2].tight_layout()
plt.show()