from funcs import *
l=1
m=101
epsilon_arr=np.linspace(0,2,m)[::-1] # Epsilon array going from 2 to 0
x_init=np.linspace(0,l,m)
y_init=np.zeros((32,m)) # Initial guess (for eps=2)
DOS_at_mid1,sol_arr1=DOS_diff_eps(epsilon_arr, x_init, y_init,l,m,phi_L=1,phi_R=0)
sol_midx=np.zeros((32,m))

for i in range(m):
    sol_midx[:,i]=sol_arr1[:,int(m/2)+i*m]

jint=current_int(sol_midx,m)
plt.plot(epsilon_arr,jint)
plt.title("Current integrand as a function of energy with $\Delta \Phi=1$")
plt.ylabel("$j(x,eps)$")
plt.xlabel("Energy, $eps$")
plt.grid()
plt.show()