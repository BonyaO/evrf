# basic_evrf.sage - Basic eVRF implementation using shared components
load("evrf_common.sage")
import time

# =============================================================================
# BASIC eVRF SPECIFIC SETUP
# =============================================================================

# Number of constraints for basic eVRF
n_basic = 4*l + 6
# Number of variables for basic eVRF
m_prime_basic = 4*l + 5
m_basic = 2^(ceil(log(m_prime_basic + n_basic, 2))) - n_basic

# Initialize generators for basic eVRF
G_vec = [find_subgroup_generator(E, q) for _ in range(n_basic + m_basic)]
H_vec = [find_subgroup_generator(E, q) for _ in range(n_basic + m_basic)]

G = find_subgroup_generator(E, q)
H = find_subgroup_generator(E, q)

G1 = G_vec[0]
GT1 = G_vec[1]
GT2 = G_vec[2]

# GT1 generator - using a fixed point for consistency
"""GT1 = E([0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb,
0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1])
assert (q*GT1).is_zero()"""

# GS generator
GS = find_subgroup_generator(F, s)

def hash_to_curve_gs(Q, x):
    """Hash function for basic eVRF. Hashes a point Q and an integer x to a point on the source group"""
    data = point_to_bytes(Q) + str(x).encode()
    # Hash the data
    h = hashlib.sha256(data).digest()

    return hash_to_curve_bandersnatch(int.from_bytes(h, byteorder='big') % q)


# Generate a random private key
def keygen():

    # Generate k according to paper
    k = ZZ.random_element(1, 2^(l + 1))

    #verification key
    Q = k * GT1

    vk = Q

    return k, vk


def generate_witness_vector(k, H_x, x_P):
    ell = l
    z = vector(Fq, n_basic)
    # define c_i constants

    c = [1] + [2] * (ell - 1)
    c.append(s - (2*ell - 1))


    # Get binary representation of k
    k_bits = Integer(k).digits(base=2, padto=ell+1)

    P = [ (2^i) * H_x for i in range(ell + 1) ]

    
    # Fill the first 3 elements of the witness vector
    z[0] = 1 
    z[1] = k

    z[2] = x_P

    delta = []
    L = []

    delta.append(k_bits[0]*H_x + c[0]*GS)
    L.append(delta[0])

    # Set k_0, x_L_0, y_L_0 in the witness
    z[3] = k_bits[0]
    z[4] = L[0][0]
    z[5] = L[0][1]

     # Compute remaining Delta_i and L_i values
    for i in range(1, ell + 1):
        # Compute Delta_i = k_i * P_i + c_i * G_S
        delta_i = k_bits[i] * P[i] + c[i] * GS
        delta.append(delta_i)
        
        # Check if Delta_i is zero (reject if so)
        if delta[i].is_zero():
            return "reject"
        
        # Compute L_i = L_(i-1) + Delta_i
        L_i = L[i-1] + delta[i]
        L.append(L_i)
        
        # Check if L_i is zero (reject if so)
        if L[i].is_zero():
            return "reject"
        
        # Check if x-coordinates are equal (reject if so)
        if L[i-1][0] == delta[i][0]:
            return "reject"
        
        # Calculate the slope s_i for the elliptic curve addition
        s_i = (L[i-1][1] - delta[i][1]) / (L[i-1][0] - delta[i][0])
        
        # Calculate base index for this element in the witness
        base_idx = 6 + (i-1)*4
        
        # Set k_i, s_i, x_L_i, y_L_i in the witness
        z[base_idx] = k_bits[i]      
        z[base_idx + 1] = s_i        
        z[base_idx + 2] = L[i][0]     
        z[base_idx + 3] = L[i][1]     
    
    # Verify final point (L_ell) matches expected final point
    assert L[ell][0] == x_P, f"Final point x-coordinate mismatch: {L[ell][0]} != {x_P}"
    
    return z

def R1CSMatrices(X):
    
    ell = l  # Use the global ell value
    
    
    # Define c_i constants
    c = [1] + [2] * (ell - 1)
    c.append(s - (2*ell - 1))
     
    # Compute Delta_i,0 and Delta_i,1 for all i
    Delta_i_0 = []
    Delta_i_1 = []
    delta_x = []
    delta_y = []
    
    for i in range(ell + 1):
        P_i = (2^i) * X
        Delta_i_0_point = c[i] * GS
        Delta_i_1_point = P_i + c[i] * GS
        
        # Check if any of the Delta points are zero
        if Delta_i_0_point.is_zero() or Delta_i_1_point.is_zero():
            return matrix(Fq, m_prime_basic, n_basic), matrix(Fq, m_prime_basic, n_basic), matrix(Fq, m_prime_basic, n_basic)
        
        # Store the coordinates
        Delta_i_0.append((Delta_i_0_point[0], Delta_i_0_point[1]))
        Delta_i_1.append((Delta_i_1_point[0], Delta_i_1_point[1]))
        
        # Calculate delta_x and delta_y
        delta_x.append(Delta_i_1_point[0] - Delta_i_0_point[0])
        delta_y.append(Delta_i_1_point[1] - Delta_i_0_point[1])
    
    # Initialize matrices
    A = matrix(Fq, m_prime_basic, n_basic)
    B = matrix(Fq, m_prime_basic, n_basic)
    C = matrix(Fq, m_prime_basic, n_basic)
    
    # Create unit vectors
    e = lambda i: vector([1 if j == i else 0 for j in range(n_basic)])
    
    # Constraint 1: k = sum_i=0^ell 2^i * k_i
    A[0] = e(0)  
    B[0] = e(1)  
    C[0] = e(3) + sum(2^i * e(6+(i-1)*4) for i in range(1, ell+1))
    
    # Constraint 2: x_P = x_L_ell
    A[1] = e(0)  
    B[1] = e(2)  
    C[1] = e(4*ell+4)  
    
    # Constraint 3: k_0 * (1 - k_0) = 0 (ensuring k_0 ∈ {0, 1})
    A[2] = e(3)  
    B[2] = e(0) - e(3)  
    C[2] = vector([0] * n_basic)  
    
    # Constraints 3+i: k_i * (1 - k_i) = 0 (ensuring k_i ∈ {0, 1})
    for i in range(1, ell+1):
        idx = 2 + i
        A[idx] = e(6+(i-1)*4) 
        B[idx] = e(0) - e(6+(i-1)*4)  
        C[idx] = vector([0] * n_basic)  
    
    # Constraint 4+ell: Verify x_L_0 is computed correctly
    A[3+ell] = e(0)  
    B[3+ell] = Delta_i_0[0][0] * e(0) + delta_x[0] * e(3)  
    C[3+ell] = e(4)  
    
    # Constraint 5+ell: Verify y_L_0 is computed correctly
    A[4+ell] = e(0)  
    B[4+ell] = Delta_i_0[0][1] * e(0) + delta_y[0] * e(3)  
    C[4+ell] = e(5)  
    
    # Constraints for each addition step
    for i in range(1, ell+1):
        # Constraint 5+ell+i: Verify s_i is computed correctly
        A[4+ell+i] = e(4*i+3) 
        B[4+ell+i] = e(4*(i-1)+4) - (Delta_i_0[i][0] * e(0) + delta_x[i] * e(6+(i-1)*4))  
        C[4+ell+i] = e(4*(i-1)+5) - (Delta_i_0[i][1] * e(0) + delta_y[i] * e(6+(i-1)*4))  
        
        # Constraint 5+2*ell+i: Verify x_L_i is computed correctly
        A[4+2*ell+i] = e(4*i+3)  
        B[4+2*ell+i] = e(4*i+3)  
        C[4+2*ell+i] = e(4*i+4) + e(4*(i-1)+4) + Delta_i_0[i][0] * e(0) + delta_x[i] * e(6+(i-1)*4)  
        
        # Constraint 5+3*ell+i: Verify y_L_i is computed correctly
        A[4+3*ell+i] = e(4*i+3)  
        B[4+3*ell+i] = e(4*(i-1)+4) - e(4*i+4)  
        C[4+3*ell+i] = e(4*(i-1)+5) + e(4*i+5)  
    
    # Pad with zeros to reach power of 2 dimension
    A_padded = matrix(Fq, m_basic, n_basic)
    B_padded = matrix(Fq, m_basic, n_basic)
    C_padded = matrix(Fq, m_basic, n_basic)
    
    # Copy the filled part
    for i in range(m_prime_basic):
        A_padded[i] = A[i]
        B_padded[i] = B[i]
        C_padded[i] = C[i]
    
    return A_padded, B_padded, C_padded


def eval(k,x):
    # Compute Q = k · G_T,1
    Q = k * GT1
    
    # Compute H_x = H(Q, x) ∈ G*_S
    H_x = hash_to_curve_gs(Q, x)
    
    # Compute P = k · H_x ∈ G*_S
    P = k * H_x
    
    # Extract x-coordinate of P
    x_P = P[0]  
    y_P = P[1]  
    
    # Set y = x_P ∈ F_q and Y = y · G_T,2
    y = x_P
    Y = y * GT2

    z = generate_witness_vector(k, H_x, x_P)

    A, B, C = R1CSMatrices(H_x)

    T = G1 + Q + Y 

    

    # Step 1: Generate Schnorr proof for Q = k*G_T1
    pi_Q = proof_R_dlog(Q, GT1, k)
    
    # Step 2: Generate Schnorr proof for Y = x_P*G_T2
    pi_Y = proof_R_dlog(Y, GT2, x_P)
    
    
    # Step 3: Generate Bulletproof for the R1CS statement
    pi_BP = generate_bulletproof(A, B, C, T, z)
    
    pi = {
        "pi_Q": pi_Q,
        "pi_Y": pi_Y,
        "pi_BP": pi_BP
    }

    return y, Y, pi


def verify(vk, x, Y, pi):
    # Decompose pi to get individual proofs
    pi_Q = pi["pi_Q"]
    pi_Y = pi["pi_Y"]
    pi_BP = pi["pi_BP"]

    # Extract Q and Y from vk
    Q = vk

    T = G1 + Q + Y
    # Step 1: Verify Schnorr proof for Q = k*G_T1 
    if not verify_proof_R_dlog(Q, GT1, pi_Q):
        return False

    # Step 2: Verify Schnorr proof for Y = x_P*G_T2
    if not verify_proof_R_dlog(Y, GT2, pi_Y):
        return False

    H_x = hash_to_curve_gs(Q, x)
    # Step 3: Verify Bulletproof for the R1CS statement
    A, B, C = R1CSMatrices(H_x)

    if not verify_bulletproof(A, B, C, T, pi_BP):
        return False


    return True


def test():
    # Generate a random private key
    k, vk = keygen()

    # Generate a random x
    x = Fq.random_element()

    # Evaluate the eVRF
    y_evrf, Y, pi = eval(k, x)

    # Verify the eVRF
    assert verify(vk, x, Y, pi), "eVRF verification failed!"
    print("eVRF verification succeeded!")


test()
