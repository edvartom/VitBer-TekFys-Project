from funcs import *
delta =0.01
eps=2

l=1 # Length of piece of metal
m=101 # Number of nodes
x=np.linspace(0,l,m) # Position
y=np.zeros((32,m)) # Initial guess

# Plotting
DOS, sol=density_of_states(fun_for_bvp, bc_for_bvp_2i, x, y, eps,l)
plt.plot(x,DOS,color="orchid")
plt.title("DOS with superconductors on either side, epsilon=2")
plt.xlabel("$x$")
plt.ylabel("$D/D_0$")
plt.grid()
plt.show()