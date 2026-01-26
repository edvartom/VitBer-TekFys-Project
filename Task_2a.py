from funcs import *
M=np.array([[1+2j,3+4j],[5+6j,7+8j]])       # Make a random complex 2x2 matrix
vector=matrix_to_vector(M)                  # Vector as specified in project description
N=vector_to_matrix(vector)                  # Transform vector back into a 2x2 matrix
result=(N==M)                               # A 2x2 matrix only containing 'True' (if equal) or 'False' (if different)
print(all(result.reshape(1,4)[0]))          # See if all elements are True
