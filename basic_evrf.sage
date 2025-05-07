# This code is a basic implementation of an Exponent Verifiable Random Function (EVRF) using the BLS12-381 curve (Target) and Bandersnatch curve (Source).
import hashlib

# base fild of BLS12-381
p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
# order of BLS12-381
q = 52435875175126190479447740508185965837690552500527637822603658699938581184513
# order bandersnatch
s = 15172417585395309745210573063711216967055694857434315578142854216712503379

assert p in Primes()
assert q in Primes()
assert s in Primes()

Fp = GF(p)
Fq = GF(q)
Fs = GF(s)

# BLS12-381 G1
E = EllipticCurve(Fp, [0, 4])
# Bandersnatch
F = EllipticCurve_from_j(Fq(8000))
a = 52435875175126190479447740508185965837690552500527637822603658699938430656513
b = 52435875175126190479447740508185965837690552500527637822603658699309173440513

# GT1 generator
GT1 = E([0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb,
0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1])
assert (q*GT1).is_zero()

# Public parameters
def calculate_ell(s, q):
    min_order = min(s, q)
    log_min_order = min_order.nbits() - 1  
    return log_min_order - 1

l = calculate_ell(E.order(), F.order())

def find_subgroup_generator(Curve, subgroup_order):
    cof1 = Curve.order() // subgroup_order
    P1 = cof1 * Curve.random_point()
    while P1.is_zero() or not (subgroup_order * P1).is_zero():
        P1 = cof1 * Curve.random_point()
    assert (subgroup_order * P1).is_zero()
    return P1

# Number of constraints
n = 4*l + 6
# Number of variables
m_prime = 4*l + 5
m = 2^(ceil(log(m_prime + n, 2))) - n

# Use find_subgroup_generator to get G_vec
G_vec = [find_subgroup_generator(E, q) for _ in range(n+m)]
# Construct H vector
H_vec = [find_subgroup_generator(E, q) for _ in range(n+m)]

G = G_vec[0]
H = H_vec[0]



# GT1 generator
GT2 = find_subgroup_generator(E, q)

# GS generator
GS = find_subgroup_generator(F, s)




def hash_to_curve_bandersnatch(x):
    while True:
        try: 
            P = F.lift_x(x)
            return P   
        except (ValueError):
            x = x+1

def point_to_bytes(P):
    """Serialize elliptic curve point to bytes in a consistent format"""
    x = Integer(P[0])
    y = Integer(P[1])
    return x.str().encode() + y.str().encode()

# Generate a random private key
def keygen():

    # Generate k according to paper
    k = ZZ.random_element(1, 2^(l + 1))

    #verification key
    Q = k * GT1

    #Schnorr proof to generate pi_Q
    pi_Q = proof_R_dlog(Q, GT1, k)

    vk = (Q, pi_Q)

    return k, vk

#Schnorr proof of knowledge of discrete log
def proof_R_dlog(P, G_T, k):

    r_nonce = Fq.random_element() #Random nonce
    R = r_nonce*G_T 

    challenge_data = point_to_bytes(P) + point_to_bytes(G_T) + point_to_bytes(R)

    challenge = Fq(int.from_bytes(hashlib.sha256(challenge_data).digest(), byteorder='big'))

    
    #Challeng response
    s = r_nonce + challenge * k

    return (R, s)


#Verify Discrete log proof that Q = k*G for some k
def verify_proof_R_dlog(P, G_T, proof):
     
    R, s = proof

    challenge_data = point_to_bytes(P) + point_to_bytes(G_T) + point_to_bytes(R)

    challenge = Fq(int.from_bytes(hashlib.sha256(challenge_data).digest(), byteorder='big'))
    
    
    left = Integer(s)*G_T
    right = R + challenge * P

    return left == right



def generate_witness_vector(k, H_x, x_P):
    ell = l
    z = vector(Fq, n)
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
            return matrix(Fq, m, n), matrix(Fq, m, n), matrix(Fq, m, n)
        
        # Store the coordinates
        Delta_i_0.append((Delta_i_0_point[0], Delta_i_0_point[1]))
        Delta_i_1.append((Delta_i_1_point[0], Delta_i_1_point[1]))
        
        # Calculate delta_x and delta_y
        delta_x.append(Delta_i_1_point[0] - Delta_i_0_point[0])
        delta_y.append(Delta_i_1_point[1] - Delta_i_0_point[1])
    
    # Initialize matrices
    A = matrix(Fq, m_prime, n)
    B = matrix(Fq, m_prime, n)
    C = matrix(Fq, m_prime, n)
    
    # Create unit vectors
    e = lambda i: vector([1 if j == i else 0 for j in range(n)])
    
    # Constraint 1: k = sum_i=0^ell 2^i * k_i
    A[0] = e(0)  # 1
    B[0] = e(1)  # k
    C[0] = e(3) + sum(2^i * e(6+(i-1)*4) for i in range(1, ell+1))
    
    # Constraint 2: x_P = x_L_ell
    A[1] = e(0)  # 1
    B[1] = e(2)  # x_P
    C[1] = e(4*ell+4)  # x_L_ell
    
    # Constraint 3: k_0 * (1 - k_0) = 0 (ensuring k_0 ∈ {0, 1})
    A[2] = e(3)  # k_0
    B[2] = e(0) - e(3)  # 1 - k_0
    C[2] = vector([0] * n)  # 0
    
    # Constraints 3+i: k_i * (1 - k_i) = 0 (ensuring k_i ∈ {0, 1})
    for i in range(1, ell+1):
        idx = 2 + i
        A[idx] = e(6+(i-1)*4)  # k_i
        B[idx] = e(0) - e(6+(i-1)*4)  # 1 - k_i
        C[idx] = vector([0] * n)  # 0
    
    # Constraint 4+ell: Verify x_L_0 is computed correctly
    A[3+ell] = e(0)  # 1
    B[3+ell] = Delta_i_0[0][0] * e(0) + delta_x[0] * e(3)  # x_Delta_0,0 + delta_x,0 * k_0
    C[3+ell] = e(4)  # x_L_0
    
    # Constraint 5+ell: Verify y_L_0 is computed correctly
    A[4+ell] = e(0)  # 1
    B[4+ell] = Delta_i_0[0][1] * e(0) + delta_y[0] * e(3)  # y_Delta_0,0 + delta_y,0 * k_0
    C[4+ell] = e(5)  # y_L_0
    
    # Constraints for each addition step
    for i in range(1, ell+1):
        # Constraint 5+ell+i: Verify s_i is computed correctly
        A[4+ell+i] = e(4*i+3)  # s_i
        B[4+ell+i] = e(4*(i-1)+4) - (Delta_i_0[i][0] * e(0) + delta_x[i] * e(6+(i-1)*4))  # x_L_(i-1) - (x_Delta_i,0 + delta_x,i * k_i)
        C[4+ell+i] = e(4*(i-1)+5) - (Delta_i_0[i][1] * e(0) + delta_y[i] * e(6+(i-1)*4))  # y_L_(i-1) - (y_Delta_i,0 + delta_y,i * k_i)
        
        # Constraint 5+2*ell+i: Verify x_L_i is computed correctly
        A[4+2*ell+i] = e(4*i+3)  # s_i
        B[4+2*ell+i] = e(4*i+3)  # s_i
        C[4+2*ell+i] = e(4*i+4) + e(4*(i-1)+4) + Delta_i_0[i][0] * e(0) + delta_x[i] * e(6+(i-1)*4)  # x_L_i + x_L_(i-1) + x_Delta_i,0 + delta_x,i * k_i
        
        # Constraint 5+3*ell+i: Verify y_L_i is computed correctly
        A[4+3*ell+i] = e(4*i+3)  # s_i
        B[4+3*ell+i] = e(4*(i-1)+4) - e(4*i+4)  # x_L_(i-1) - x_L_i
        C[4+3*ell+i] = e(4*(i-1)+5) + e(4*i+5)  # y_L_(i-1) + y_L_i
    
    # Pad with zeros to reach power of 2 dimension
    A_padded = matrix(Fq, m, n)
    B_padded = matrix(Fq, m, n)
    C_padded = matrix(Fq, m, n)
    
    # Copy the filled part
    for i in range(m_prime):
        A_padded[i] = A[i]
        B_padded[i] = B[i]
        C_padded[i] = C[i]
    
    return A_padded, B_padded, C_padded

def hash_to_curve_gs(Q, x):
   
    data = point_to_bytes(Q) + str(x).encode()
    # Hash the data
    h = hashlib.sha256(data).digest()
    
    # Convert the hash to an integer
    h_int = int.from_bytes(h, byteorder='big')

    return hash_to_curve_bandersnatch(h_int % q)
    

# Helper function to get Fiat-Shamir challenges
def compute_challenge(data):
    """Compute a challenge using the Fiat-Shamir heuristic"""
    if isinstance(data, list):
        data = b''.join([str(d).encode() for d in data])
    elif not isinstance(data, bytes):
        data = str(data).encode()
    return Fq(int.from_bytes(hashlib.sha256(data).digest(), byteorder='big'))

def debug_inner_product(u, v, location):
    """Helper function to debug inner product calculations"""
    result = inner_product(u, v)
    # Append inner product debug info to a file
    with open('inner_product_debug.txt', 'a') as f:
        f.write(f"\nInner product at {location}:\n")
        f.write(f"Vector u : {[u[i] for i in range(len(u))]}\n")
        f.write(f"Vector v : {[v[i] for i in range(len(v))]}\n") 
        f.write(f"Result: {result}\n")
    return result

def print_vector(v):
    with open('vector_debug.txt', 'a') as f:
        f.write(f"\nVector:\n")
        f.write(f"{[v[i] for i in range(len(v))]}\n")

    return 


def generate_bulletproof(A, B, C, T, z):

    m_b, n_b = A.nrows(), A.ncols()
    ell = l  # Using the global l value from the original code

    # In this implementation, z' is a zero vector and eta is 0 in Fq
    z_prime = vector(Fq, n_b, [0] * n_b)
    eta = Fq(0)

    # r is the size of x in the witness vector z = (x||y)
    r = 3
    x = vector(Fq, [z[i] for i in range(r)])
    y = vector(Fq, [z[i] for i in range(r, n_b)])
    
    # Compute Az
    Az = A * z
    
    # Compute Bz
    Bz = B * z

    
    # Create zero vector of length n_b for the first part of the H term
    zero_n = vector(Fq, [0] * n_b)
    
    # Compute the first part: ⟨((x′||y) || Az), G⟩
    # Since x' is all zeros (z_prime is all zeros), we only need (0r||y)
    x_prime = vector(Fq, [0] * r)  # This is x′ (all zeros)
    xy_concat = vector(Fq, list(x_prime) + list(y))  # Concatenate x′ and y
    G_term_vec = vector(Fq, list(xy_concat) + list(Az))  # Concatenate (x′||y) with Az
    G_term = 0
    for i in range(len(G_term_vec)):
        G_term += G_term_vec[i] * G_vec[i]
    
    # Compute the second part: ⟨(0n|| Bz), H⟩
    H_term_vec = vector(Fq, list(zero_n) + list(Bz))  # Concatenate 0n with Bz
    H_term = 0
    for i in range(len(H_term_vec)):
        H_term += H_term_vec[i] * H_vec[i]
    
    # Compute the final commitment S
    S = G_term + H_term + r * H  # H_base is the single generator H from the paper
    
    # For Fiat-Shamir, serialize S
    S_bytes = point_to_bytes(S)

    # Step 3: Use Fiat-Shamir to get the verifier's challenges
    transcript_so_far = [s]
    alpha = compute_challenge(transcript_so_far + [b"alpha"])
    beta = compute_challenge(transcript_so_far + [b"beta"])
    gamma = compute_challenge(transcript_so_far + [b"gamma"])
    delta = compute_challenge(transcript_so_far + [b"delta"])


    # Update transcript
    transcript_so_far.extend([alpha, beta, gamma, delta])
    
    # Step 4: Calculate intermediate values as per the paper
    mu = alpha * gamma

   

    # calculate delta vector as in the protocol
    delta_vec = vector(Fq, n_b)
    for i in range(r):
        delta_vec[i] = delta
    for i in range(r, n_b):
        delta_vec[i] = 1

    
    # Calculate delta_inv vector
    delta_inv_vec = vector(Fq, n_b)
    for i in range(r):
        delta_inv_vec[i] = delta^(-1)
    for i in range(r, n_b):
        delta_inv_vec[i] = 1

    
    

    # Construct G' vector with modified generators
    G_prime = []
    for i in range(n_b):
        G_prime.append(G_vec[i])
    for i in range(m_b):
        G_prime.append(gamma^(-1) * G_vec[n_b+i])
    

    # Compute vector c from the matrices
    
    mu_vec = vector(Fq, [mu for i in range(m_b)])
    beta_vec = vector(Fq, [beta for i in range(m_b)])
    gamma_vec = vector(Fq, [gamma for i in range(m_b)])

    


    A_mu = mu_vec * A
    B_beta = beta_vec * B
    C_gamma = gamma_vec * C

    
    c = A_mu + B_beta - C_gamma

    
    
    # Compute inner product value omega
    omega = 0
    # Create vectors for first term
    alpha_vec = vector(Fq, [alpha for i in range(m_b)])
    beta_vec = vector(Fq, [beta for i in range(m_b)])
    omega = inner_product(alpha_vec, beta_vec)

    # Create vectors for second term 
    alpha_vec2 = vector(Fq, [alpha for i in range(n_b)])
    h_product = hadamard_product(c, delta_vec)

    

    delta_c_product = inner_product(alpha_vec2, h_product)
 




    omega += delta^2 * delta_c_product

    
    # Compute the public point P
    P = delta_inv_vec[0] * T + S
    
    # Add the G' term
    G_prime_term = 0
    for i in range(n_b):
        G_prime_term += delta^2 * alpha^(i+1) * G_prime[i]
    for i in range(m_b):
        G_prime_term -= beta^(i+1) * G_prime[n_b+i]
    P += G_prime_term
    
    # Add the H term
    H_term = 0
    for i in range(n_b):
        H_term += c[i] * delta_vec[i] * H_vec[i]
    for i in range(m_b):
        H_term -= alpha^(i+1) * H_vec[n_b+i]
    P += H_term

     # Step 5: The prover computes vectors u and v for the inner product argument
    # Initialize u vector
    u = vector(Fq, n_b+m_b)
    
    # First part (x'||y) + delta^(-1) * (x||y') + delta^2 * alpha^n_b
    for i in range(r):
        u[i] = z_prime[i] + delta_inv_vec[i] * z[i] + delta^2 * alpha^(i+1)
    for i in range(r, n_b):
        u[i] = z[i] + delta_inv_vec[i] * z_prime[i] + delta^2 * alpha^(i+1)
    
    # Second part (Az + delta^(-1) * Az') * gamma^m_b - beta^m_b
    Az_prime = A * z_prime  # This is a zero vector since z_prime is zero
    for i in range(m_b):
        u[n_b+i] = (Az[i] + delta_inv_vec[0] * Az_prime[i]) * gamma^(i+1) - beta^(i+1)
    
    # Initialize v vector
    v = vector(Fq, n_b+m_b)
    
    # First part c * delta
    for i in range(n_b):
        v[i] = c[i] * delta_vec[i]
    
    # Second part Bz - alpha^m_b + delta^(-1) * Bz'
    Bz_prime = B * z_prime  # This is a zero vector since z_prime is zero
    for i in range(m_b):
        v[n_b+i] = Bz[i] - alpha^(i+1) + delta_inv_vec[0] * Bz_prime[i]
    
    # Compute eta'
    eta_prime = r + delta_inv_vec[0] * eta
    
    # Initialize the inner product argument
    inner_product_proof = generate_inner_product_proof(G_prime, H_vec, G, H, P, omega, u, v, eta_prime)
    
    # Return the complete bulletproof
    return {
        "S": S,
        "challenges": {
            "alpha": alpha,
            "beta": beta,
            "gamma": gamma,
            "delta": delta
        },
        "inner_product_proof": inner_product_proof
    }
def generate_inner_product_proof(G, H, G_base, H_base, P, omega, a, b, alpha):

    
    challenge_data = point_to_bytes(G_base) + point_to_bytes(H_base) + point_to_bytes(P) + str(omega).encode()
    e = Fq(int.from_bytes(hashlib.sha256(challenge_data).digest(), byteorder='big'))
   
    
    
    G_prime = e * G_base
    P_Prime = P + omega * G_prime

    
    return generate_modified_inner_product_proof(G, H, G_prime, H_base, P_Prime, a, b, alpha)


def generate_modified_inner_product_proof(G, H, G_base, H_base, P, a, b, alpha):
    n = len(a)
    if n == 1:
        # For n=1, we simply reveal the values a, b, and alpha
        r = Fq.random_element()
        s = Fq.random_element()
        delta = Fq.random_element()
        eta = Fq.random_element()
        
        # Compute commitments A and B.
        #TODO: Confirm with Dr. Atonio if this is correct according to the paper
        A = r * G_base + s * H_base + (r * b[0] + s * a[0]) * G_base + delta * H_base
        B = (r * s) * G_base + eta * H_base
        
        # Generate Fiat-Shamir challenge
        e = compute_challenge([A, B])
        
        # Compute responses
        r_prime = r + a[0] * e
        s_prime = s + b[0] * e
        delta_prime = eta + delta * e + alpha * e^2
        
        return {
            "A": A,
            "B": B,
            "e": e,
            "r_prime": r_prime,
            "s_prime": s_prime,
            "delta_prime": delta_prime
        }
    
    # Recursive case
    n_half = n // 2
    
    # Split vectors
    a_L = a[:n_half]
    a_R = a[n_half:]
    b_L = b[:n_half]
    b_R = b[n_half:]
    G_L = G[:n_half]
    G_R = G[n_half:]
    H_L = H[:n_half]
    H_R = H[n_half:]
    
    # Compute inner products for the left and right sides
    c_L = inner_product(a_L, b_R)
    c_R = inner_product(a_R, b_L)
    
    # Sample random blinding factors
    d_L = Fq.random_element()
    d_R = Fq.random_element()
    
    # Compute commitments L and R
    L = inner_product(a_L, G_R) + inner_product(b_R, H_L) + c_L * G_base + d_L * H_base
    R = inner_product(a_R, G_L) + inner_product(b_L, H_R) + c_R * G_base + d_R * H_base
    
    # Generate Fiat-Shamir challenge
    x = compute_challenge([L, R])

    x_inv = pow(x, -1, q)

    # Compute new generators and vectors for the next level
    G_prime = [x_inv * g + x * h for g, h in zip(G_L, G_R)]
    H_prime = [x * g + x_inv * h for g, h in zip(H_L, H_R)]
    a_prime = [x * a + x_inv * b for a, b in zip(a_L, a_R)]
    b_prime = [x_inv * a + x * b for a, b in zip(b_L, b_R)]
    
  
    # Compute the new point P
    P_prime = x^2 * L + P + pow(x, -2, q) * R
    
    # Compute new alpha
    alpha_prime = x^2 * d_L + alpha + pow(x, -2, q) * d_R
    
    # Recursively generate the inner product proof
    inner_proof = generate_modified_inner_product_proof(G_prime, H_prime, G_base, H_base, P_prime, a_prime, b_prime, alpha_prime)
    
    # Combine the current level with the recursive proof
    return {
        "L": L,
        "R": R,
        "x": x,
        "inner_proof": inner_proof
    }
def inner_product(u, v):
    if len(u) != len(v):
        raise ValueError("Vectors must have the same length")
    
    if hasattr(u, 'parent') and hasattr(v, 'parent'):
        # For Sage vectors, use built-in dot product
        return u.dot_product(v)
    else:
        # For regular sequences, use manual computation
        result = 0
        for i in range(len(u)):
            result += u[i] * v[i]
        return result

def hadamard_product(u, v):
    if len(u) != len(v):
        raise ValueError("Vectors must have the same length")
    
    if hasattr(u, 'parent') and hasattr(v, 'parent'):
        # Check if they're Sage vectors over the same field
        return vector(u.parent().base_ring(), [u[i] * v[i] for i in range(len(u))])
    else:
        # Generic implementation for lists or other sequence types
        return [u[i] * v[i] for i in range(len(u))]


def verify_bulletproof(A, B, C, T, proof):

    # Extract proof components
    S = proof["S"]
    challenges = proof["challenges"]
    alpha = challenges["alpha"]
    beta = challenges["beta"]
    gamma = challenges["gamma"]
    delta = challenges["delta"]
    inner_product_proof = proof["inner_product_proof"]

    m, n = A.nrows(), A.ncols()
    ell = l  # Using the global l value from the original code
    
    r = 3 # r is the size of x in the witness vector z = (x||y)

    
    # Calculate intermediate values as in the proof generation
    mu = alpha * gamma

    # Create delta vector
    delta_vec = vector(Fq, n)
    for i in range(r):  # r should be defined or passed as parameter
        delta_vec[i] = delta
    for i in range(r, n):
        delta_vec[i] = 1
    
    
    # Calculate delta_inv vector
    delta_inv_vec = vector(Fq, n)
    for i in range(r):
        delta_inv_vec[i] = delta^(-1)
    for i in range(r, n):
        delta_inv_vec[i] = 1

    # Construct G' vector with modified generators
    G_prime = []
    for i in range(n):
        G_prime.append(G_vec[i])
    for i in range(m):
        G_prime.append(gamma^(-1) * G_vec[n+i])
    
    
    # Compute vector c from the matrices
    mu_vec = vector(Fq, [mu for i in range(m)])
    beta_vec = vector(Fq, [beta for i in range(m)])
    gamma_vec = vector(Fq, [gamma for i in range(m)])

    
    
    A_mu = mu_vec * A
    B_beta = beta_vec * B
    C_gamma = gamma_vec * C


    c = A_mu + B_beta - C_gamma

   
    
    # Compute inner product value omega
    omega = 0
    # Create vectors for first term
    alpha_vec = vector(Fq, [alpha for i in range(m)])
    beta_vec = vector(Fq, [beta for i in range(m)])
    omega = inner_product(alpha_vec, beta_vec)

    # Create vectors for second term 
    alpha_vec2 = vector(Fq, [alpha for i in range(n)])
    h_product = hadamard_product(c, delta_vec)

    
    
    delta_c_product = inner_product(alpha_vec2, h_product)


    omega += delta^2 * delta_c_product

    # Compute the public point P
    P = delta_inv_vec[0] * T + S
    
    # Add the G' term
    G_prime_term = 0
    for i in range(n):
        G_prime_term += delta^2 * alpha^(i+1) * G_prime[i]
    for i in range(m):
        G_prime_term -= beta^(i+1) * G_prime[n+i]
    P += G_prime_term
    
    # Add the H term
    H_term = 0
    for i in range(n):
        H_term += c[i] * delta_vec[i] * H_vec[i]
    for i in range(m):
        H_term -= alpha^(i+1) * H_vec[n+i]
    P += H_term

    # Verify the inner product proof
    return verify_inner_product_proof(G_prime, H_vec, G, H, P, omega, inner_product_proof)

def verify_inner_product_proof(G, H, G_base, H_base, P, omega, proof):

   
    challenge_data = point_to_bytes(G_base) + point_to_bytes(H_base) + point_to_bytes(P) + str(omega).encode()
    e = Fq(int.from_bytes(hashlib.sha256(challenge_data).digest(), byteorder='big'))

    G_prime = e * G_base

    
    P_Prime = P + omega * G_prime
    

    return verify_modified_inner_product_proof(G, H, G_prime, H_base, P_Prime, proof)

def verify_modified_inner_product_proof(G, H, G_base, H_base, P, proof):
    # Base case: vectors of length 1
    n = len(G)
    if n == 1:
        # Extract the components of the final proof
        A = proof["A"]
        B = proof["B"]
        e = proof["e"]
        r_prime = proof["r_prime"]
        s_prime = proof["s_prime"]
        delta_prime = proof["delta_prime"]
        
        # The correct verification equation should be:
        # e^2 * P + e * A + B = (r_prime * e * G_base) + (s_prime * e * H_base) + (r_prime * s_prime * G_base) + (delta_prime * H_base)
        
        left_side = e^2 * P + e * A + B
        
        # Build right side term by term
        term1 = r_prime * e * G_base      # r'eG
        term2 = s_prime * e * H_base      # s'eH
        term3 = r_prime * s_prime * G_base # r's'G
        term4 = delta_prime * H_base       # δ'H
        
        right_side = term1 + term2 + term3 + term4
        

        return left_side == right_side
    
    # Recursive case
    # Extract components from the proof
    L = proof["L"]
    R = proof["R"]
    x = proof["x"]
    inner_proof = proof["inner_proof"]
    
    # Verify the challenge matches what we'd compute
    if x != compute_challenge([L, R]):
        return False
    
    # Calculate the new values for the recursive verification
    n_half = n // 2
    G_L = G[:n_half]
    G_R = G[n_half:]
    H_L = H[:n_half]
    H_R = H[n_half:]
    
    x_inv = pow(x, -1, q)

    # Compute new generators and vectors for the next level
    G_prime = [x_inv * g + x * h for g, h in zip(G_L, G_R)]
    H_prime = [x * g + x_inv * h for g, h in zip(H_L, H_R)]
    
    # Compute the new point P'
    P_prime = x^2 * L + P + pow(x, -2, q) * R
    
    # Recursively verify the inner product proof
    return verify_modified_inner_product_proof(G_prime, H_prime, G_base, H_base, P_prime, inner_proof)


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

    T = GT1 + Q + Y 

    # Step 1: Generate Schnorr proof for Q = k*G_T1
    pi_Q = proof_R_dlog(Q, GT1, k)
    
    # Step 2: Generate Schnorr proof for Y = x_P*G_T2
    pi_Y = proof_R_dlog(Y, GT2, x_P)
    
    # Step 3: Generate Bulletproof for the R1CS statement
    pi_BP = generate_bulletproof(A, B, C, T, z)
    
    # Combine all
   
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
    Q = vk[0]
    Y_vk = vk[1]

    #Define commitment T
    T = GT1 + Q + Y
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

    # Compute the public key Y
    y, Y, pi = eval(k, x)

    # Verify the proof
    assert verify(vk, x, Y, pi), "Verification failed"
    print("Verification successful")


    # Compute the hash to curve for H_x
    H_x = hash_to_curve_bandersnatch(x)

test()
