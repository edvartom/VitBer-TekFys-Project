from funcs import *     # Secant method function is in funcs (the original, didn't really change it)

# Defining quantities for RK32
h0 = 0.005
tol = 1e-7
alpha = 0.8
x_init = 0              # Starting position
x_end = 2*np.pi         # End position
b0 = 0 # Guess 1
b1 = 8 # Guess 2

prev_b=[]               # List for storing b-values that have already been plotted.
                        # Comes in handy for the function below.

# Making a function with RK32 so that it suits the format of the secant method   
def RK32_for_secant(b):
    '''
    INPUT
        b: Float. y'(0)=b. Secant_method must find a b so that y(2*pi)=0
    
    OUTPUT
        y_arr[-1,0]: Float. Last y-value calculated by RK32. Equivalent to y(2*pi).
                     When this is close enough to zero the secant method is finished.
    '''
    y_init = np.array([0,b],dtype=np.float64)  # y(x_init)=0, y'(x_init)=b
    x_arr, y_arr, h_arr, N = RK32(x_init, x_end, y_init, f_1c, h0, tol, alpha) # Calculate y with RK32
    if b not in prev_b:  # If y(x) with this b has not been plotted yet...
        plt.plot(x_arr,y_arr[:,0], label='b ='+str(b))  # ...then plot it...
        prev_b.append(b)                                # ... and store b-value.
    return y_arr[-1,0]                                  # =y(2*pi).

b = secant_method(b0, b1, RK32_for_secant, tol=0.01) # Calculate b with secant method
upper_bound = RK32_for_secant(b) # =y(2*pi)

#Plotting 
plt.title("Calculations of boundary value problem y(2pi)=0 and y'(0)=b")
plt.xlabel('$x$')
plt.ylabel('$y(x)$')
plt.plot(2*np.pi,0, 'go', label='Boundary cond.') # Marking where we want end point to be
plt.grid()
plt.legend()
plt.show()
print("b = ",b)
print('y(2*pi) = ',upper_bound)

