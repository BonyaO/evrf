# Full DDH eVRF implementation using shared components
load("evrf_common.sage")

# =============================================================================
# Full DDH eVRF SPECIFIC SETUP
# =============================================================================


# Number of constraints
n = 7*l + 8
# Number of variables
m_prime = 7*l + 7
m = 2^(ceil(log(m_prime + n, 2))) - n

# Use find_subgroup_generator to get G_vec
G_vec = [find_subgroup_generator(E, q) for _ in range(n+m)]
# Construct H vector
H_vec = [find_subgroup_generator(E, q) for _ in range(n+m)]

G = find_subgroup_generator(E, q)
H = find_subgroup_generator(E, q)

G1 = G_vec[0]
GT1 = G_vec[1]
GT2 = G_vec[2]

# GS generator
GS = find_subgroup_generator(F, s)

# Hash function
def hash_to_curve_gs(Q, x, j):
   
    data = point_to_bytes(Q) + str(x).encode() + str(j).encode()
    # Hash the data
    h = hashlib.sha256(data).digest()
    
    # Convert the hash to an integer
    h_int = int.from_bytes(h, byteorder='big')

    return hash_to_curve_bandersnatch(h_int % q)


# Generate a random private key
def keygen():

    k = ZZ.random_element(1, 2^(l + 1))

    # Sample a random k' in Fq
    k_prime = Fq.random_element()
    
    # Verification key
    Q = k * GT1

    vk = (Q, k_prime)
    sk = (k, k_prime)

    return sk, vk


def generate_witness_vector(y, X1, X2, k):
    # Global parameters
    ell = l  
    z = vector(Fq, n)

    # define c_i constants
    c = [1] + [2] * (ell - 1)
    c.append(s - (2*ell - 1))
    
    # Compute binary representation of k
    k_bits = Integer(k).digits(base=2, padto=ell+1)
    
    # Initialize z with first 3 elements (statement part)
    z[:3] = [1, k, y]

    G_S = GS

    # Fill in k_bits (binary representation of k)
    for i in range(len(k_bits)):
        z[3 + i] = k_bits[i]

    # Starting index for witness data after k_bits
    witness_idx = 3 + len(k_bits)  # This should be 3 + (ell + 1)
    
    for j in [1, 2]:  # Process both X1 and X2
        X = X1 if j == 1 else X2
        
        # Compute P_i = 2^i · X for i in {0, ..., ℓ}
        P_i_values = []
        Delta_i_values = []
        P_i_0 = X 
        P_i_values.append(P_i_0)
        Delta_0 = k_bits[0] * P_i_0 + c[0] * G_S
        Delta_i_values.append(Delta_0)
        
        # Check if Δ_0 = 0, if so reject
        if Delta_0.is_zero():
            return "reject"
        
        # Compute L_i points
        L_i_values = []
        L_i = Delta_i_values[0]  
        L_i_values.append(L_i)
        
        # Check if L_0 = 0, if so reject
        if L_i.is_zero():
            return "reject"
        
        # Add x_L_0, y_L_0 to witness
        z[witness_idx] = L_i[0]      
        z[witness_idx + 1] = L_i[1]  
        witness_idx += 2
        
        s_i_values = []

        for i in range(1, ell + 1):
            P_i = (2^i) * X
            P_i_values.append(P_i)

            Delta_i = k_bits[i] * P_i_values[i] + c[i] * G_S
            Delta_i_values.append(Delta_i)
            
            # Check if Δ_i = 0, if so reject
            if Delta_i.is_zero():
                return "reject"

            L_i = L_i_values[i-1] + Delta_i_values[i]
            L_i_values.append(L_i)
            
            # Check if L_i = 0, if so reject
            if L_i.is_zero():
                return "reject"
            
            # Check if x-coordinates are equal, if so reject
            if L_i_values[i-1][0] == Delta_i_values[i][0]:
                return "reject"
            
            # Compute s_i values
            x_L_prev = L_i_values[i-1][0]
            y_L_prev = L_i_values[i-1][1]
            x_Delta_i = Delta_i_values[i][0]
            y_Delta_i = Delta_i_values[i][1]

            s_i = (y_L_prev - y_Delta_i) * ((x_L_prev - x_Delta_i)^(-1)) % q
            s_i_values.append(s_i)
            
            # Add w_i^(j) = (s_i^(j), x_L_i^(j), y_L_i^(j)) to witness
            z[witness_idx] = s_i        
            z[witness_idx + 1] = L_i[0] 
            z[witness_idx + 2] = L_i[1] 
            witness_idx += 3
       
    return z

def R1CSMatricesFull(X1, X2, k_prime):
    # Global parameters
    ell = l  # This should be defined globally as ⌊log₂ min{s, q}⌋ - 1
    G_S = GS  # Generator of G_S
    
    # Define c_i constants
    c = [1] + [2] * (ell - 1)
    c.append(s - (2*ell - 1))
    
    # Precompute Delta values for both X^(1) and X^(2)
    Delta_values = {}
    
    for j in [1, 2]:
        X = X1 if j == 1 else X2
        Delta_values[j] = {}
        
        for i in range(ell + 1):
            # Compute P_i^(j) = 2^i * X^(j)
            P_i = (2^i) * X
            
            # Compute Delta_{i,0}^(j) = c_i * G_S
            Delta_i_0 = c[i] * G_S
            
            # Compute Delta_{i,1}^(j) = P_i + c_i * G_S
            Delta_i_1 = P_i + c[i] * G_S
            
            # Check if any Delta points are zero
            if Delta_i_0.is_zero() or Delta_i_1.is_zero():
                # Return zero matrices as specified in the paper
                A = matrix(Fq, m, n)
                B = matrix(Fq, m, n) 
                C = matrix(Fq, m, n)
                return (A, B, C)
            
            # Extract coordinates
            x_delta_0, y_delta_0 = Delta_i_0[0], Delta_i_0[1]
            x_delta_1, y_delta_1 = Delta_i_1[0], Delta_i_1[1]
            
            # Store Delta values
            Delta_values[j][i] = {
                'x_0': x_delta_0,
                'y_0': y_delta_0,
                'x_1': x_delta_1,
                'y_1': y_delta_1,
                'delta_x': x_delta_1 - x_delta_0,
                'delta_y': y_delta_1 - y_delta_0
            }
    
    # Initialize matrices
    A = matrix(Fq, m_prime, n)
    B = matrix(Fq, m_prime, n)
    C = matrix(Fq, m_prime, n)
    
    # Unit vector helper
    e = lambda i: vector([1 if j == i else 0 for j in range(n)])
    
    # Build the constraint matrices row by row
    
    # Row 1: k = sum_{i=0}^ell 2^i * k_i
    A[0] = e(0)  
    B[0] = e(1)  
    C[0] = sum(2^i * e(3+i) for i in range(ell+1))
    
    # Row 2: y = k' * x_L_ell^(1) + x_L_ell^(2)
    A[1] = e(0)  
    B[1] = e(2)  
    C[1] = k_prime * e(4+4*ell) + e(6+7*ell)
    
    # Rows 3+i: k_i * (1 - k_i) = 0 for i in {0, ..., ell}
    for i in range(ell + 1):
        A[2+i] = e(3+i)  
        B[2+i] = e(0) - e(3+i)  
        C[2+i] = vector([0] * n)  
    
    # Current row index
    row = 3 + ell
    
    # Constraints for j=1 and j=2
    for j in [1, 2]:
        # Row 3+ell+j: x_L_0^(j) constraint
        A[row] = e(0)  
        B[row] = Delta_values[j][0]['x_0'] * e(0) + Delta_values[j][0]['delta_x'] * e(3)
        C[row] = e(4+ell+(j-1)*(3*ell+2))  
        row += 1
    
    # Row 5+ell+j: y_L_0^(j) constraint
    for j in [1, 2]:
        A[row] = e(0)  
        B[row] = Delta_values[j][0]['y_0'] * e(0) + Delta_values[j][0]['delta_y'] * e(3)
        C[row] = e(5+ell+(j-1)*(3*ell+2))  
        row += 1
    
    # Constraints for each iteration i in [1, ..., ell]
    for j in [1, 2]:
        for i in range(1, ell + 1):
            # Row 5+j+2i+ell: s_i^(j) constraint
            A[row] = e(3+ell+3*i+(j-1)*(3*ell+2)) 
            B[row] = (e(4+ell+3*(i-1)+(j-1)*(3*ell+2)) - 
                     Delta_values[j][i]['x_0'] * e(0) - 
                     Delta_values[j][i]['delta_x'] * e(3+i))
            C[row] = (e(5+ell+3*(i-1)+(j-1)*(3*ell+2)) - 
                     Delta_values[j][i]['y_0'] * e(0) - 
                     Delta_values[j][i]['delta_y'] * e(3+i))
            row += 1
    
    # x_L_i constraints
    for j in [1, 2]:
        for i in range(1, ell + 1):
            # Row 5+j+2i+3*ell: x_L_i^(j) constraint
            A[row] = e(3+ell+3*i+(j-1)*(3*ell+2))  
            B[row] = e(3+ell+3*i+(j-1)*(3*ell+2))  
            C[row] = (e(4+ell+3*i+(j-1)*(3*ell+2)) + 
                     e(4+ell+3*(i-1)+(j-1)*(3*ell+2)) + 
                     Delta_values[j][i]['x_0'] * e(0) + 
                     Delta_values[j][i]['delta_x'] * e(3+i))
            row += 1
    
    # y_L_i constraints
    for j in [1, 2]:
        for i in range(1, ell + 1):
            # Row 5+j+2i+5*ell: y_L_i^(j) constraint
            A[row] = e(3+ell+3*i+(j-1)*(3*ell+2))  
            B[row] = (e(4+ell+3*(i-1)+(j-1)*(3*ell+2)) - 
                     e(4+ell+3*i+(j-1)*(3*ell+2)))
            C[row] = (e(5+ell+3*(i-1)+(j-1)*(3*ell+2)) + 
                     e(5+ell+3*i+(j-1)*(3*ell+2)))
            row += 1
    
    # Pad matrices to reach the required dimensions
    A_padded = matrix(Fq, m, n)
    B_padded = matrix(Fq, m, n)
    C_padded = matrix(Fq, m, n)
    
    # Copy the filled part
    for i in range(m_prime):
        A_padded[i] = A[i]
        B_padded[i] = B[i] 
        C_padded[i] = C[i]
    
    return (A_padded, B_padded, C_padded)


def eval(sk,x):

    k, k_prime = sk
    
    # Compute Q = k * GT1
    Q = k * GT1
    
    # Compute the two hash points
    X1 = hash_to_curve_gs(Q, x, 1)
    X2 = hash_to_curve_gs(Q, x, 2)
    
    # Compute P1 and P2
    P1 = k * X1
    P2 = k * X2
    
    # Extract x-coordinates
    x_P1 = P1[0]
    x_P2 = P2[0]
    
    # Compute the final output y using the leftover hash lemma extraction
    y = (k_prime * x_P1 + x_P2) % q
    
    Y = y * GT2
    
    # First, generate proofs for the discrete logarithm relations
    pi_Q = proof_R_dlog(Q, GT1, k)
    pi_Y = proof_R_dlog(Y, GT2, y)
    
    
    z = generate_witness_vector(y, X1, X2, k)

    # Generate R1CS matrices
    A, B, C = R1CSMatricesFull(X1, X2, k_prime)


    # Generate bulletproof for both R1CS statements
    pi_BP = generate_bulletproof(A, B, C, T, z)
    
    # Combine all proofs
    pi = {
        "pi_Q": pi_Q,
        "pi_Y": pi_Y,
        "pi_BP1": pi_BP,
    }
    return y, Y, pi


def verify(vk, x, Y, pi):
    """Verify an output of the full DDH-based eVRF"""
    Q, k_prime = vk
    pi_Q = pi["pi_Q"]
    pi_Y = pi["pi_Y"]
    pi_BP1 = pi["pi_BP1"]
    
    # Hash the input to get the challenge points
    X1 = hash_to_curve_gs(Q, x, 1)
    X2 = hash_to_curve_gs(Q, x, 2)


    # Verify the discrete logarithm proofs
    if not verify_proof_R_dlog(Q, GT1, pi_Q):
        return False

    if not verify_proof_R_dlog(Y, GT2, pi_Y):
        return False

    # Verify the bulletproof for the R1CS statement
    A, B, C = R1CSMatricesFull(X1, X2, k_prime)
    
    T = G1 + Q + Y
    
    return verify_bulletproof(A, B, C, T, pi_BP1)
    
    
# Test the implementation
def test():
    print("Testing keygen...")
    sk, vk = keygen()
    k, k_prime = sk
    Q, k_prime_pub = vk

    x = Fq.random_element()

    # Evaluate VRF
    y, Y, pi = eval(sk, x)
    
    # Verify the output
    assert verify(vk, x, Y, pi), "VRF verification failed"
    print("Eval and verify test passed!")


test()