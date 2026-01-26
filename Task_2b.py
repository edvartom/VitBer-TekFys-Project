from funcs import *

# Defining vectors
m1 = np.array([1,2,3,4,1,2,3,4])
m2 = np.array([5,6,7,8,5,6,7,8])
m3 = np.array([9,10,11,12,9,10,11,12])
m4 = np.array([13,14,15,16,13,14,15,16])

# Transforming
v = m_vectors_to_v(m1, m2, m3, m4)
print(v)
m_reconstructed = v_to_m_vectors(v)
print(m_reconstructed)

# Checking and printing result
assert np.allclose((m1, m2, m3, m4), m_reconstructed), "m_vectors_to_v and v_to_m_vectors do not execute inverse transformations"
print("m_vectors_to_v and v_to_m_vectors are indeed inverse transformations")