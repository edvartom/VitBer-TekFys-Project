from funcs import *

l = 1                       # Metal piece length
m = 101                     # Number of nodes
x = np.linspace(0,l,m)      # Position
y = np.zeros((32, m))       # Initial guess
ε_arr = np.array([0, 1, 2]) # Different ε guesses
tol = 0.0001                # Tolerance

for eps in ε_arr:  # Finding solution wih solve_bvp for each ε
    # λ function allows us to change ε and l in fun and bc
    sol = solve_bvp(lambda x, vec: fun_for_bvp(x, vec, ε = eps), 
                    bc_for_bvp_2f, x, y)
    not_zero = 0  # Checking how many of the elements in y are non-zero.

    # Set a tolerance limit. With numerical inaccuracy it may not be EXACTLY zero 
    for i in range(32):
        for j in range(m):
            if abs(sol.y[i, j]) > tol:  # Iterate through every element in y and 
                                        # checking whether the absolute value is
                                        # greater than the tolerance
                not_zero += 1

    print("ε = ", eps)
    print("Number of elements in solution that are not numerically equal to zero: ", not_zero,"\n")
    # not_zero = 0, as expected for normal metals