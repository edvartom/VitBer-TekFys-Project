from funcs import *
m=101
l=1
epsilon_arr=np.linspace(0,2,m)[::-1] # Epsilon array going from 2 to 0
eps_of_int=np.array([2.0,1.5,1.0,0.5,0.0])
x_arr=np.linspace(0,l,m)
y_init=np.zeros((32,m)) # Initial guess (for eps=2)
DOS_at_mid1,sol_arr1=DOS_diff_eps(epsilon_arr, x_arr, y_init,l,m) #This can be deleted when in jupyter notebooks

plt.title("Current integrand $j(x,eps)$ for $eps={2.0,1.5,1.0,0.5,0.0}$")
for epsilon in eps_of_int:
    ind=int(np.where(epsilon_arr==epsilon)[0])
    jint=current_int(sol_arr1[:,m*ind:(m+m*ind)],m) # This is from 2k
    plt.plot(x_arr,jint,label="epsilon="+str(epsilon))
plt.ylabel("$j(x,eps)$")
plt.xlabel("Position, $x$")
plt.grid()
plt.legend()
plt.show()