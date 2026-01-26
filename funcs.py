# We shall define all functions here for now
# Later We will put this in jupyter file. 
# Use "import funcs as *" to use the functions directly.

import numpy as np
import numpy.typing as npt
from collections.abc import Callable
from scipy.integrate import solve_bvp
from scipy.integrate import simpson
import matplotlib.pyplot as plt

#Task 1c
def y_anl_1(x: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    '''
    Analytic solution for vector y in task 1.
    INPUT
        x: (?,) NDArray (1D array with x values)
    OUTPUT
        y: (?,2) NDArray (2D array with the analytic solution [y(x), y'(x)])
    '''
    return np.transpose(np.array([np.sin(2*x), 2*np.cos(2*x)]))

def f_1c(x: float, y: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    '''
    OUTPUT
        y: (2,) NDArray (1D array)
    '''
    return np.array([y[1], -4*np.sin(2*x)])

# Defining functions

def y_anl_1(x: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    '''
    Analytic solution for the vector y in task 1.
    INPUT
        x:      (?,) NDArray (1D array with x values)
    OUTPUT
        y:      (?,2) NDArray (2D array with the analytic 
                solution [y(x), y'(x)])
    '''
    return np.transpose(np.array([np.sin(2*x), 2*np.cos(2*x)]))

def f_1c(x: float, y: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    '''
    OUTPUT
        y:  (2,) NDArray (1D array)
    '''
    return np.array([y[1], -4*np.sin(2*x)])

def RK32(x_init: float, x_end: float, y_init: npt.NDArray[np.float64], 
         f: Callable, h0: float, tol: float, α: float
         ) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], 
        npt.NDArray[np.float64], npt.NDArray[np.float64], int]: 
    '''
    INPUT
        y_init:     2D array of initial values for y: [y(x_init), y'(x_init)]
        f:          the derivative of y
        h0:         initial step length
        tol:        maximum error tolerated
        α:      pessimist factor to make the step lenght h a bit 
                    smaller than calculated
    OUTPUT
        y: (?,?) NDArray (2D array with y(x_n) on row n)
        x: (?,?) NDArray (1D array with x-values)
        h: (?,) NDArray (1D array of step lengths)
        N: int (total numbers of steps taken)
    ''' 
    h, N = h0, 0
    x_lst, y_lst, h_lst = [x_init], [y_init], [h0]
    k1 = f(x_init, y_init)

    while (x_end - x_lst[-1]) > 1e-15:                              # In other words: while x_n < x_end
        N += 1                                                      # Count number of loops
        x, y = x_lst[-1], y_lst[-1].copy()                          # We use the last list element for calculating the subsequent
        h = np.min([h, x_end-x])                                    # Ensures we do not calculate y-values for x>x_endh=np.min(h,x_end-x)                             # Assures we do not calculate y-values for x>x_end
        k2 = f(x + 0.5*h, y + 0.5* h * k1)
        k3 = f(x + 0.75*h, y + 0.75* h * k2)
        y += (1/9) * h * (2*k1 + 3*k2 + 4*k3)                       # Update y
        x += h                                                      # Update x
        k4 = f(x, y)                                                # k4 uses the new values for x and y
        z = y_lst[-1] + (1/24) * h*(7*k1 + 6*k2 + 8*k3 + 3*k4)
        est = np.linalg.norm(y-z)                                   # Calculate the error (norm)
        if est<tol:                                                 # If the error is tolerated, we may use the values for y, (and update k1)
            k1 = k4 
            x_lst.append(x)
            y_lst.append(y.copy()) 
            h_lst.append(h)                                         # k4 can be reused as k1 in the next loop
        h *= α*(tol/est)**(1/3)                                     # h is updated with the same formula whether or not the values were accepted
    return np.array(x_lst), np.array(y_lst), np.array(h_lst), N     # Returning numpy arrays

#Task 1d
def error_afo_tol(x_init: float, x_end: float, y_init: npt.NDArray[np.float64], 
                  f: Callable, h0: float, tol: npt.NDArray[np.float64], 
                  α: float, f_anl: Callable) -> npt.NDArray[np.float64]:
    '''
    Calculates the max error for a goiven function y

    INPUT
        y_init:     2D array of initial values for y: [y(x_init), y'(x_init)]
        f:          the derivative of y
        h0:         initial step length
        tol:        array of maximum error tolerated
        α:      pessimist factor to make the step length h a bit smaller 
                    than calculated
        y_anl:      analytic function for y
    OUTPUT
        err_arr:      Maximum error array for different tolerances
    '''
    err_arr = np.zeros(len(tol))
    for i in range(len(tol)):
        x_arr, y_arr, h_arr, N = RK32(x_init, x_end, y_init, f, h0, tol[i], α)
        y_anl = f_anl(x_arr)
        err_max = np.max(np.abs((y_arr[:,0] - y_anl[:,0])))
        err_arr[i] = err_max
    return err_arr

def N_timesteps_afo_α(x_init: float, x_end: float, 
                          y_init: npt.NDArray[np.float64], f: Callable, 
                          h0: float, tol: float, α: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    '''
    INPUT
        y_init:     2D array of initial values for y: [y(x_init), y'(x_init)]
        f:          the derivative of y
        h0:         initial step length
        tol:        maximum error tolerated
        α:      array of pessimist factors to make the step lenght h a 
                    bit smaller than calculated
    OUTPUT
        N: time steps for each α
    '''
    N_timesteps_array = np.zeros(len(α))
    for i in range(len(α)):
        x_arr, y_arr, h_arr, N_timesteps = RK32(x_init, x_end, y_init, f, h0, tol, α[i]) # The function count the number of times the while loop is running
        N_timesteps_array[i] = (N_timesteps)
    return N_timesteps_array

# Task 1e
def secant_method(z0: float, z1: float, g: Callable, tol=0.001):
    '''
    Finds root of function g

    INPUT
        z0: float, initial guess 1
        z1: float, initial guess 2
        g: function. The function for which we want to find a root
        tol: float, tolerance
    OUTPUT
        z[-1]: float. The root of g, the z value for which g(z[-1]) = 0
    '''
    z = [z0, z1] 

    while tol < np.abs(z[-2] - z[-1]):
        if g(z[-1]) - g(z[-2]) == 0:  # Avoid division by zero
            print("Secant method failed (division by zero).")
            return None
        
        zi = (z[-2] * g(z[-1]) - z[-1] * g(z[-2])) / (g(z[-1]) - g(z[-2]))
        z.append(zi)

    return z[-1]

# Task 1g

def f_1g(x: float, y: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    '''
    OUTPUT
        y: (2,) NDArray (1D array)
    '''
    return np.array([y[1], y[0] + np.sin(x)])

def BVP_solver(root_finder, IVP_solver, b0, b1, x_init, x_end, f, h0, tol, α):
    '''
    Solves BVPs as an IVP, given y_init=0

    INPUT:
        root_finder:    Callable that finds x for which y(x)=0
        IVP_solver:     Callable that solves IVP
        b0:             initial guess for root  
        b1:             second initial guess for root
        f:              function that derivates y for IVP
        h0:             initial step length in IVP
        tol:            maximum error tolerated in IVP
        α:          pessimist factor to adjust step length in IVP

    OUTPUT: 
        y: (?,?) NDArray (2D array with y(x_n) on row n)
        x: (?,?) NDArray (1D array with x-values)
        h: (?,) NDArray (1D array of step lengths)
        N: int (total number of steps taken in final solution)
    '''
    
    def IVP_solver_root_finder(b):

        y_init = np.array([0,b],dtype=np.float64)
        x_arr, y_arr, h_arr, N = IVP_solver(x_init, x_end, y_init, 
                                            f, h0, tol, α)

        return y_arr[-1,0]
        
    b = root_finder(b0, b1, IVP_solver_root_finder, tol=0.01)
    y_init_new = np.array([0,b],dtype=np.float64)
    
    return IVP_solver(x_init, x_end, y_init_new, f, h0, tol, α)

#Task 2a
def matrix_to_vector(matrix):
    '''
    Transforms a 2x2 complex matrix into a 8-component real vector as 
    specified in task 2a.

    INPUT
        matrix: (2,2) N_array
    
    OUTPUT
        vector: (1,8) N_array
    
    '''
    assert matrix.shape == (2,2)
    vector = np.zeros(8)
    flat_matrix =matrix.reshape(1,4)[0]  # Reshaping matrix into (1,4) array 
                                         # with ax=0
    vector[0:4] = np.real(flat_matrix)
    vector[4:8] = np.imag(flat_matrix)
    return vector

def vector_to_matrix(vector):
    '''
    Does the exact same thing as the unction above only backwards. 
    Transforms one 8-component vector into a (2x2) matrix
    '''
    assert vector.size == 8
    matrix=np.zeros((2,2), dtype=np.complex_)   # Default dtype is float. 
                                                # Convert to complex
    matrix[0,0] = complex(vector[0],vector[4])  # a=complex(1,2) -> a=1+2.j
    matrix[0,1] = complex(vector[1],vector[5])
    matrix[1,0] = complex(vector[2],vector[6])
    matrix[1,1] = complex(vector[3],vector[7])
    return matrix

#Task 2b
def m_vectors_to_v(m1, m2, m3, m4):
    '''
    Turns four 8-component vectors m1, m2, m3, m4 into one 32-component vector v.
    Assumption: The four first numbers in the m-vectors are real values from 
    matrix in task 2a.
                When making v, the real values are grouped in the first 16 slots, 
                and the imaginary are grouped in the last 16 slots.
    Test code:
        m1 = np.array([1, 2, 3, 4, 1, 2, 3, 4])
        m2 = np.array([5, 6, 7, 8, 5, 6, 7, 8])
        m3 = np.array([9, 10, 11, 12, 9, 10, 11, 12])
        m4 = np.array([13, 14, 15, 16, 13, 14, 15, 16])
        print( m_vectors_to_v(m1, m2, m3, m4) )
    
    '''
    assert (m1.size, m2.size, m3.size, m4.size) == (8, 8, 8, 8) # Ensuring that all arrays have size 8
    v = np.concatenate((m1[0:4], m2[0:4], m3[0:4], m4[0:4],
                        m1[4:8], m2[4:8], m3[4:8], m4[4:8]))   # Merging arrays as explained above
    return v

def v_to_m_vectors(v):
    '''
    Does the exact opposite of the function above. Turns a 32-component vector v into 
    four 8-component vectors m
    Test code:
        v=np.array(range(0,32))
        print(v_to_m_vectors(v))
    
    '''
    assert v.size == 32
    m1 = np.concatenate((v[0:4], v[16:20]))
    m2 = np.concatenate((v[4:8], v[20:24]))
    m3 = np.concatenate((v[8:12], v[24:28]))
    m4 = np.concatenate((v[12:16], v[28:32]))
    return m1, m2, m3, m4

#Task 2c
def matrices_to_v(ɣ, ɣ_tilde, ω, ω_tilde):
    '''
    Converts four (2x2) matrices into a 32-component vector
    '''
    assert (ɣ.shape, ɣ_tilde.shape, ω.shape, 
            ω_tilde.shape) == ((2, 2),(2, 2),(2, 2),(2, 2))
    
    g1 = matrix_to_vector(ɣ)
    g2 = matrix_to_vector(ɣ_tilde)
    o1 = matrix_to_vector(ω)
    o2 = matrix_to_vector(ω_tilde)

    v = m_vectors_to_v(g1, g2, o1, o2)
    return v

def v_to_matrices(v):
    '''
    Converts a 32-component vector into four (2x2) matrices into four matrixes:
    (a_r, b_r, c_r, d_r, e_r, f_r, g_r, h_r, j_r, k_r, l_r, m_r, n_r, o_r, p_r, q_r,
     a_i, b_i, c_i, d_i, e_i, f_i, g_i, h_i, j_i, k_i, l_i, m_i, n_i, o_i, p_i, q_i)
     
     ---->

     ((a_r+i*a_i, b_r+i*b_i), (c_r+i*c_i, d_r+i*d_i)),
     ((e_r+i*e_i, f_r+i*f_i), (g_r+i*g_i, h_r+i*h_i)),
     ((j_r+i*j_i, k_r+i*k_i), (l_r+i*l_i, m_r+i*m_i)),
     ((n_r+i*n_i, o_r+i*o_i), (p_r+i*p_i, q_r+i*p_i))
    '''
    assert v.size == 32
    g1, g2, o1, o2 = v_to_m_vectors(v)
    ɣ = vector_to_matrix(g1)
    ɣ_tilde = vector_to_matrix(g2)
    ω = vector_to_matrix(o1)
    ω_tilde = vector_to_matrix(o2)

    return ɣ, ɣ_tilde, ω, ω_tilde


#Task 2d
δ = 0.01 # Constant 

def N_matrix(ɣ, ɣ_tilde):
    '''
    Creates the matrix N as specified in project description.

    INPUT
        ɣ:       (2,2) N_array
        ɣ_tilde: (2,2) N_array
    
    OUTPUT
        N:           (2,2) N_array
    '''
    # Making all matrices complex to avoid problems with dtype
    ɣ, ɣ_tilde = ɣ.astype(complex), ɣ_tilde.astype(complex)

    I = np.array([[1,0], [0,1]]).astype(complex) # Identity matrix
    N = I - ɣ @ ɣ_tilde
    a, b, c, d = N[0,0], N[0,1], N[1,0], N[1,1]
    N = ( 1 / (a*d - b*c) ) * np.array([[d, -b], [-c, a]]) #Inverting the matrix
    return N

def N_tilde_matrix(ɣ, ɣ_tilde):
    '''
    Creates the matrix N_tilde as specified in project description.
    This is exactly the same function as N_matrix except that ɣ and ɣ_tilde 
    switch places in einsum
    '''
    ɣ, ɣ_tilde = ɣ.astype(complex), ɣ_tilde.astype(complex)

    I = np.array([[1,0],[0,1]]).astype(complex)
    N = I - ɣ_tilde @ ɣ
    a, b, c, d = N[0,0], N[0,1], N[1,0], N[1,1]
    N = ( 1 / (a*d - b*c) ) * np.array([[d, -b],[-c, a]]) #Inverting the matrix
    return N

def diff_v(v,ε):
    '''
    Differentiates vector v. 
    Converts v into four matrices (ɣ, ɣ_tilde, ω, ω_tilde) and these are 
    differentiated as per eq. (9)-(12) in the project description. The 
    differentiated matrices are then transformed back into a 32-component vector dv.
    Test code:
        ε = 2
        ɣ = np.array([[1+2j, 3+4j], [5+6j, 7+8j]])
        ɣ_tilde = np.array([[9+10j, 11+12j], [13+14j, 15+16j]])
        ω = np.array([[17+18j, 19+20j], [21+22j, 23+24j]])
        ω_tilde = np.array([[25+26j,27+28j],[29+30j,31+32j]])
        v = matrices_to_v(ɣ,ɣ_tilde,ω,ω_tilde)
        dv = diff_v(v,ε)
        print(dv)


    INPUT (do not change input variables, they are specified in project desc.)
        v:        (1,32) N_array
        ε:  Float. Dimensionless energy quantity

    OUTPUT
        dv:      (1,32) N_array. Derivative of v
    
    '''
    ɣ, ɣ_tilde, ω, ω_tilde = v_to_matrices(v)
    N, N_tilde = N_matrix(ɣ, ɣ_tilde), N_tilde_matrix(ɣ, ɣ_tilde)
    dω = complex(0,-2) * complex(ε,δ) * ɣ - 2 * ω @ N_tilde @ ɣ_tilde @ ω
    dω_tilde = complex(0,-2) * complex(ε, δ) * ɣ_tilde - 2 * ω_tilde @ N @ ɣ @ω_tilde
    dv = matrices_to_v(ω, ω_tilde, dω, dω_tilde)
    return dv

#Task 2e
def fun_for_bvp(x, vec, ε = 2):
    '''
    INPUT
        x:      (1,m) NDArray. Positions in interval [0,l] where l is the length 
                of the metal that we are studying.

        vec:    (32,m) NDArray. m different v-vectors for each position x.
                A vector v represents four "rolled out" matrices: 
                ɣ, ɣ_tilde, ω, ω_tilde, which is what we are trying to find.

    O UTPUT
        dv_arr: (32,m) NDArray. An array of the derivatives of the v-vectors in vec.
    '''
    
    dv_arr = np.zeros((vec.shape[0], vec.shape[1]))
    for i in range (vec.shape[1]):
        dv_arr[:,i] = diff_v(vec[:,i], ε)  # Differentiate each column (vector v)
    return dv_arr


#Task 2f

def bc_for_bvp_2f(v_left, v_right, ε = 2, length = 1, φ_L =0, φ_R=0):
    '''
    Calculates how far off the solution is from the boundary conditions 
    (eq. 13-16 in project desc.)

    INPUT
        v_left:  (1,32) NDArray. Riccati matrices (ɣ) and their derivatives 
        (ω) in x=0
        v_right: (1,32) NDArray. Riccati matrices (ɣ) and their derivatives 
        (ω) in x=l
    
    OUTPUT
        boundary_condition_res: (1,32) NDArray. Residuals of boundary conditions.
                                Ideally the entire array should be zero, 
                                but numerical solutions may not be perfect.
    
    '''
    ζ = 3 # Constant, as per project description

    ɣ_L = np.zeros((2,2))
    ɣ_tilde_L = np.zeros((2,2))
    ɣ_R = np.zeros((2,2))
    ɣ_tilde_R = np.zeros((2,2))

    #Folding up v_left and v_right into matrices
    ɣ_0, ɣ_tilde_0, ω_0, ω_tilde_0 = v_to_matrices(v_left)
    ɣ_l, ɣ_tilde_l, ω_l, ω_tilde_l = v_to_matrices(v_right)

    I = np.array([[1,0], [0,1]]).astype(complex) # Identity matrix

    #Calculating N-matrices
    N_L = N_matrix(ɣ_L, ɣ_tilde_L)
    N_tilde_L = N_tilde_matrix(ɣ_L,ɣ_tilde_L)
    N_R = N_matrix(ɣ_R, ɣ_tilde_R)
    N_tilde_R = N_tilde_matrix(ɣ_R, ɣ_tilde_R)

    #Partsums to avoid long lines:
    psum_ω_0 = (1 / (ζ * length)) * (I - ɣ_0 @ ɣ_tilde_L)
    psum_ω_l = (1 / (ζ * length)) * (I - ɣ_l @ ɣ_tilde_R)
    psum_ω_tilde_0 = (1 / (ζ * length)) * (I - ɣ_tilde_0 @ ɣ_L)
    psum_ω_tilde_l = (1 / (ζ * length)) * (I - ɣ_tilde_l @ ɣ_R)

    #Calculating eq (13)-(16)
    eq_13 = ω_0 + psum_ω_0 @ N_L @ (ɣ_L - ɣ_0)
    eq_14 = ω_tilde_0 + psum_ω_tilde_0 @ N_tilde_L @ (ɣ_tilde_L - ɣ_tilde_0)
    eq_15 = ω_l - psum_ω_l @ N_R @ (ɣ_R - ɣ_l)
    eq_16 = ω_tilde_l - psum_ω_tilde_l @ N_tilde_R @ (ɣ_tilde_R - ɣ_tilde_l)

    boundary_condition_res = matrices_to_v(eq_13, eq_14, eq_15, eq_16)
    return boundary_condition_res

# Task 2h
def density_of_states(fun, bc, x, y, eps,l,φ_L=0,φ_R=0):
    '''
    Calcuates the density of states (DOS) as per eq. (19) in project description.

    INPUT
        fun:     Callable. Function for solve_bvp
        bc:      Callable. Evaluation of boundary conditions for solve_bvp
        x:       (1,m) NDArray. Positions/nodes
        y:       (32,m) NDArray. Initial guess
        ε: Double. Dimensionless energy constant

    OUTPUT
        DOS:     (1,32) NDArray. Density of states for different x
        sol:     (32,m) NDArray. Solution calculated by bvp
    '''
    # Finding solution wih solve_bvp
    # λ function allows us to change ε and l in fun and bc
    sol = solve_bvp(lambda x, vec: fun(x,vec,ε=eps), lambda v_left, 
                    v_right: bc(v_left, v_right, ε = eps, length = l,
                                φ_L = φ_L, φ_R = φ_R), x, y)
    
    nodes = sol.x   # Positions
    m = nodes.size  # Number of nodes
    vecs = sol.y    # m vectors of size 32

    I = np.array([[1, 0],[0, 1]]).astype(complex)       # Identity matrix
    zero = np.array([[0, 0],[0, 0]]).astype(complex)    # Zero matrix
    ρ3 = np.block([[I, zero],[zero, -I]])   # 4x4 matrix, specified in project desc. 
                                            # np.block makes one matrix out of several

    DOS = np.zeros(m)                 # DOS = Density of states (value for every x)
    for i in range (m):               # Iterating through every vector in solution
        ɣ, ɣ_tilde, ω, ω_tilde = v_to_matrices(vecs[:,i]) # Unpacking each vector into four matrices
        
        # Defining N matrices:
        N = N_matrix(ɣ,ɣ_tilde)
        N_tilde = N_tilde_matrix(ɣ,ɣ_tilde)

        # Defining each element of Green function and putting them together into one matrix:
        green00 = 2 * N - I
        green01 = 2 * N @ ɣ
        green10 = -2 * N_tilde @ ɣ_tilde
        green11 = -2 * N_tilde + I
        Green = np.block([[green00, green01], [green10, green11]])

        # Calculating DOS as per eq. (19):
        DOS[i] = np.real(np.trace(ρ3 @ Green)) / 4
    return DOS, sol

#Task 2i
δ = 0.01 # Constant 

# Defining functions from project description:
def θ_sig(σ, ε):
    '''
    INPUT
        σ:   Int. 1 or -1
        ε: Double. Dimensionless energy constant
    
    OUTPUT
        Complex
    '''
    return np.arctanh(σ / complex(ε, δ))

def s_sig(σ, ε):
    '''
    INPUT
        σ:   Int. 1 or -1
        ε: Double. Dimensionless energy constant
    OUTPUT
        Complex
    '''
    return np.sinh(θ_sig(σ, ε))

def c_sig(σ, ε):
    '''
    INPUT
        σ:   Int. 1 or -1
        ε: Double. Dimensionless energy constant
    OUTPUT
        Complex
    '''
    return np.cosh(θ_sig(σ,ε))

φ_L, φ_R = 0,0 # Phase of the left and right superconductor

# Defining Riccati functions of superconductors as per project description:
def interface_Riccati(ε,φ,φ_R):
    '''
    INPUT
        ε: Double. Dimensionless energy constant
        φ_L:   Double. Phase of left superconductor
        φ_R:   Double. Phase of right superconductor

    OUTPUT
        ɣ_L:        (2,2) NDArray. Riccati function of left superconductor
        ɣ_tilde_L:  (2,2) NDArray. Riccati function of left superconductor
        ɣ_R:        (2,2) NDArray. Riccati function of right superconductor
        ɣ_tilde_L:  (2,2) NDArray. Riccati function of right superconductor

    '''
    s_plus = s_sig(1, ε) / (1 + c_sig(1, ε))
    s_minus = s_sig(-1, ε) / (1 + c_sig(-1, ε))

    ɣ_L = np.array([[0, s_plus], [s_minus, 0]]) * np.exp(complex(0, φ_L))
    ɣ_tilde_L = np.array([[0, s_minus], [s_plus, 0]]) * np.exp(complex(0, -φ_L))
    ɣ_R = np.array([[0, s_plus], [s_minus, 0]]) * np.exp(complex(0, φ_R))
    ɣ_tilde_R = np.array([[0, s_minus], [s_plus, 0]]) * np.exp(complex(0, -φ_R))
    
    return ɣ_L, ɣ_tilde_L, ɣ_R, ɣ_tilde_R


# making new bc for superconductor interfaces
def bc_for_bvp_2i(v_left,v_right,ε=2,length=1,phi_L=0,phi_R=0):
    '''
    Calculates how far off the solution is from the boundary conditions (eq. 13-16 in project desc.)

    INPUT
        v_left:  (1,32) NDArray. Riccati matrices (ɣ) and their derivatives (ω) in x=0
        v_right: (1,32) NDArray. Riccati matrices (ɣ) and their derivatives (ω) in x=l
    
    OUTPUT
        boundary_condition_res: (1,32) NDArray. Residuals of boundary conditions.
                                Ideally the entire array should be zero, 
                                but numerical solutions may not be perfect.
    
    '''
    ζ = 3 #Constant, as per project desc.

    ɣ_L, ɣ_tilde_L, ɣ_R, ɣ_tilde_R = interface_Riccati(ε, φ_L, φ_R)

    # Folding up v_left and v_right into matrices
    ɣ_0, ɣ_tilde_0, ω_0, ω_tilde_0 = v_to_matrices(v_left)
    ɣ_l, ɣ_tilde_l, ω_l, ω_tilde_l = v_to_matrices(v_right)

    I = np.array([[1,0], [0,1]]).astype(complex) # Identity matrix

    #Calculating N-matrices
    N_L = N_matrix(ɣ_L, ɣ_tilde_L)
    N_tilde_L = N_tilde_matrix(ɣ_L, ɣ_tilde_L)
    N_R = N_matrix(ɣ_R, ɣ_tilde_R)
    N_tilde_R = N_tilde_matrix(ɣ_R, ɣ_tilde_R)

    #Partsums to avoid long lines:
    psum_ω_0 = (1 / (ζ * length)) * (I - ɣ_0 @ ɣ_tilde_L)
    psum_ω_l = (1 / (ζ * length)) * (I - ɣ_l @ ɣ_tilde_R)
    psum_ω_tilde_0 = (1 / (ζ * length)) * (I - ɣ_tilde_0 @ ɣ_L)
    psum_ω_tilde_l = (1 / (ζ * length)) * (I - ɣ_tilde_l @ ɣ_R)

    #Calculating eq (13)-(16)
    eq_13 = ω_0 + psum_ω_0 @ N_L @ (ɣ_L - ɣ_0)
    eq_14 = ω_tilde_0 + psum_ω_tilde_0 @ N_tilde_L @ (ɣ_tilde_L - ɣ_tilde_0)
    eq_15 = ω_l-psum_ω_l @ N_R @ (ɣ_R - ɣ_l)
    eq_16 = ω_tilde_l - psum_ω_tilde_l @ N_tilde_R @ (ɣ_tilde_R - ɣ_tilde_l)

    boundary_condition_res=matrices_to_v(eq_13, eq_14, eq_15, eq_16)
    return boundary_condition_res

#Task 2k
def DOS_diff_eps(ε_arr, x_init, y_init,l,m,phi_L=0,phi_R=0):
    DOS_at_mid=np.zeros(ε_arr.size)
    sol_arr=np.zeros((32,m*ε_arr.size))
    for i in range (ε_arr.size): # Iterating through all εs
        DOS, sol=density_of_states(fun_for_bvp, bc_for_bvp_2i, x_init, y_init, ε_arr[i],l,phi_L=phi_L,phi_R=phi_R) #Calculating DOS for each ε
        DOS_at_mid[i]=DOS[int(x_init.size/2)] #Take value in the midddle
        vecs=sol.y[:,0:m]
        y_init=vecs.copy() # Solution is new inital guess for the next ε
        sol_arr[:,i*m:(m+i*m)]=vecs
    return DOS_at_mid, sol_arr

#Task 2l
def current_int(sol,m): 
    jint_arr=np.zeros(m)

    I=np.array([[1,0],[0,1]]).astype(complex) # Identity matrix
    Zero=np.array([[0,0],[0,0]]).astype(complex) # Zero matrix
    rho_3=np.block([[I,Zero],[Zero,-I]]) # 4x4 matrix, specified in project desc.

    for i in range (m):                  # Iterating through every vector in solution
        ɣ, ɣ_tilde, ω, ω_tilde=v_to_matrices(sol[:,i]) # Unpacking each vector into four matrices
        
        # Defining N matrices:
        N=N_matrix(ɣ,ɣ_tilde)
        dxN=N@(ω@ɣ_tilde+ɣ@ω_tilde)@N

        N_tilde=N_tilde_matrix(ɣ,ɣ_tilde)
        dxN_tilde=N_tilde@(ω_tilde@ɣ+ɣ_tilde@ω)@N_tilde

        # Defining each element of Green function and putting them together into one matrix:
        green00=2*N-I
        green01=2*N@ɣ
        green10=-2*N_tilde@ɣ_tilde
        green11=-2*N_tilde+I
        Green=np.block([[green00,green01],[green10,green11]])

        dxgreen01=N@ω+dxN@ɣ
        dxgreen10=-N_tilde@ω_tilde-dxN_tilde@ɣ_tilde
        dxGreen=np.block([[dxN,dxgreen01],[dxgreen10,-dxN_tilde]])
        jint=np.real(np.trace(rho_3@(Green@dxGreen-dxGreen@Green)))
        jint_arr[i]=jint
    return jint_arr
