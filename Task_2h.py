from funcs import*

l = 1                       # Length of piece of metal
m = 101                     # Number of nodes
x = np.linspace(0, l, m)    # Position
y = np.zeros((32, m))       # Initial guess
ε_arr = np.array([0, 1, 2]) # Try for different εs 


# Plotting
for ε in ε_arr:
    DOS, sol=density_of_states(fun_for_bvp, bc_for_bvp_2f, x, y, ε,l)
    plt.plot(x,DOS, label='ε = '+str(ε))
    plt.xlabel('Positions, $x$')
    plt.ylabel('$D(x,eps)/D_0$')
    plt.title('Density of states as a function of position for ε$={0,1,2}$')
    plt.grid()
    plt.legend()

plt.show()