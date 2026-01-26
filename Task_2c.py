from funcs import *

# Defining four complex matrices
gamma = np.array([[1+2j,3+4j],[5+6j,7+8j]])
gamma_tilde = np.array([[9+10j,11+12j],[13+14j,15+16j]])
omega = np.array([[17+18j,19+20j],[21+22j,23+24j]])
omega_tilde = np.array([[25+26j,27+28j],[29+30j,31+32j]])

# Applying transformation funtions
v = matrices_to_v(gamma,gamma_tilde,omega,omega_tilde)

M1,M2,M3,M4 = v_to_matrices(v)

# Check if inverse, print result.
assert np.allclose(M1,M2,M3,M4, v_to_matrices(v)), "Functions are not inverse transformations"
print("The functions are inverse transformations.")
print(f"32-component vector after transormation: {v}") # Writes down the components of vector v