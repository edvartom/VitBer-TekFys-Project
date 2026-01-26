from funcs import *

# phi_R should be 0 all the time
phi_L_arr=np.linspace(0,2*np.pi,10)
l=1
m=101
epsilon_arr_2n=np.linspace(0,2,m)[::-1] # Epsilon array going from 2 to 0
x_init_2n=np.linspace(0,l,m)
y_init_2n=np.zeros((32,m)) # Initial guess (for eps=2)

j_arr = np.zeros(len(phi_L_arr))
for i in range(len(phi_L_arr)):
    DOS_at_mid1_2n,sol_arr1_2n=DOS_diff_eps(epsilon_arr_2n, x_init_2n, y_init_2n,l,m,phi_L_arr[i],phi_R=0)
    sol_midx_2n=sol_arr1_2n[:,int(m/2)::m]
    jint_2n=current_int(sol_midx_2n,m)
    j_arr[i] = - simpson(jint_2n)

fig, ax_2n = plt.subplots()
ax_2n.plot(phi_L_arr,j_arr)
ax_2n.set_title("Current as a function of phase difference")
ax_2n.set_ylabel("Currents, $I/I_0$")
ax_2n.set_xlabel("Angle difference, $\Delta \phi$")
ax_2n.grid()
plt.show()