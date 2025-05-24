# evrf_common.sage - Shared components for both basic and full eVRF implementations
import hashlib

# =============================================================================
# GLOBAL PARAMETERS AND CURVES
# =============================================================================

# Base field of BLS12-381
p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
# Order of BLS12-381
q = 52435875175126190479447740508185965837690552500527637822603658699938581184513
# Order bandersnatch
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

# Calculate ell parameter
def calculate_ell(s, q):
    min_order = min(s, q)
    log_min_order = min_order.nbits() - 1  
    return log_min_order - 1

l = calculate_ell(E.order(), F.order())

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def find_subgroup_generator(Curve, subgroup_order):
    """Find a generator of a subgroup of given order"""
    cof1 = Curve.order() // subgroup_order
    P1 = cof1 * Curve.random_point()
    while P1.is_zero() or not (subgroup_order * P1).is_zero():
        P1 = cof1 * Curve.random_point()
    assert (subgroup_order * P1).is_zero()
    return P1

def point_to_bytes(P):
    """Serialize elliptic curve point to bytes in a consistent format"""
    x = Integer(P[0])
    y = Integer(P[1])
    return x.str().encode() + y.str().encode()

def hash_to_curve_bandersnatch(x):
    """Hash to a point on the Bandersnatch curve"""
    while True:
        try: 
            P = F.lift_x(x)
            return P   
        except (ValueError):
            x = x+1

def compute_challenge(data):
    """Compute a challenge using the Fiat-Shamir heuristic"""
    if isinstance(data, list):
        data = b''.join([str(d).encode() for d in data])
    elif not isinstance(data, bytes):
        data = str(data).encode()
    return Fq(int.from_bytes(hashlib.sha256(data).digest(), byteorder='big'))

# =============================================================================
# SCHNORR PROOF SYSTEM
# =============================================================================

def proof_R_dlog(P, G_T, k):
    """Generate Schnorr proof of knowledge of discrete log"""
    r_nonce = Fq.random_element()
    R = r_nonce * G_T 

    challenge_data = point_to_bytes(P) + point_to_bytes(G_T) + point_to_bytes(R)
    challenge = Fq(int.from_bytes(hashlib.sha256(challenge_data).digest(), byteorder='big'))
    
    # Challenge response
    s = r_nonce + challenge * k
    return (R, s)

def verify_proof_R_dlog(P, G_T, proof):
    """Verify Schnorr proof that Q = k*G for some k"""
    R, s = proof

    challenge_data = point_to_bytes(P) + point_to_bytes(G_T) + point_to_bytes(R)
    challenge = Fq(int.from_bytes(hashlib.sha256(challenge_data).digest(), byteorder='big'))
    
    left = Integer(s) * G_T
    right = R + challenge * P

    return left == right

# =============================================================================
# INNER PRODUCT UTILITIES
# =============================================================================

def inner_product(u, v):
    """Compute inner product of two vectors"""
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
    """Compute Hadamard (element-wise) product of two vectors"""
    if len(u) != len(v):
        raise ValueError("Vectors must have the same length")
    
    if hasattr(u, 'parent') and hasattr(v, 'parent'):
        # Check if they're Sage vectors over the same field
        return vector(u.parent().base_ring(), [u[i] * v[i] for i in range(len(u))])
    else:
        # Generic implementation for lists or other sequence types
        return [u[i] * v[i] for i in range(len(u))]

# =============================================================================
# BULLETPROOF SYSTEM (SHARED COMPONENTS)
# =============================================================================

def generate_inner_product_proof(G, H, G_base, H_base, P, omega, a, b, alpha):
    """Generate inner product proof"""
    challenge_data = point_to_bytes(G_base) + point_to_bytes(H_base) + point_to_bytes(P) + str(omega).encode()
    e = Fq(int.from_bytes(hashlib.sha256(challenge_data).digest(), byteorder='big'))
    
    G_prime = e * G_base
    P_Prime = P + omega * G_prime
    
    return generate_modified_inner_product_proof(G, H, G_prime, H_base, P_Prime, a, b, alpha)

def generate_modified_inner_product_proof(G, H, G_base, H_base, P, a, b, alpha):
    """Generate modified inner product proof (recursive)"""
    n = len(a)
    if n == 1:
        # Base case: vectors of length 1
        r = Fq.random_element()
        s = Fq.random_element()
        delta = Fq.random_element()
        eta = Fq.random_element()
        
        # Compute commitments A and B
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

def verify_inner_product_proof(G, H, G_base, H_base, P, omega, proof):
    """Verify inner product proof"""
    challenge_data = point_to_bytes(G_base) + point_to_bytes(H_base) + point_to_bytes(P) + str(omega).encode()
    e = Fq(int.from_bytes(hashlib.sha256(challenge_data).digest(), byteorder='big'))

    G_prime = e * G_base
    P_Prime = P + omega * G_prime
    
    return verify_modified_inner_product_proof(G, H, G_prime, H_base, P_Prime, proof)

def verify_modified_inner_product_proof(G, H, G_base, H_base, P, proof):
    """Verify modified inner product proof (recursive)"""
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

def generate_bulletproof(A, B, C, T, z):
    """Generate bulletproof for R1CS statement (shared implementation)"""
    m_b, n_b = A.nrows(), A.ncols()

    # In this implementation, z' is a zero vector and eta is 0 in Fq
    z_prime = vector(Fq, n_b, [0] * n_b)
    eta = Fq(0)

    # r is the size of x in the witness vector z = (x||y)
    r = 3
    x = vector(Fq, [z[i] for i in range(r)])
    y = vector(Fq, [z[i] for i in range(r, n_b)])
    
    # Compute Az and Bz
    Az = A * z
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
    
    # Use Fiat-Shamir to get the verifier's challenges
    transcript_so_far = [s]
    alpha = compute_challenge(transcript_so_far + [b"alpha"])
    beta = compute_challenge(transcript_so_far + [b"beta"])
    gamma = compute_challenge(transcript_so_far + [b"gamma"])
    delta = compute_challenge(transcript_so_far + [b"delta"])

    # Update transcript
    transcript_so_far.extend([alpha, beta, gamma, delta])
    
    # Calculate intermediate values as per the paper
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

    # The prover computes vectors u and v for the inner product argument
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

def verify_bulletproof(A, B, C, T, proof):
    """Verify bulletproof for R1CS statement (shared implementation)"""
    # Extract proof components
    S = proof["S"]
    challenges = proof["challenges"]
    alpha = challenges["alpha"]
    beta = challenges["beta"]
    gamma = challenges["gamma"]
    delta = challenges["delta"]
    inner_product_proof = proof["inner_product_proof"]

    m, n = A.nrows(), A.ncols()
    
    r = 3 # r is the size of x in the witness vector z = (x||y)
    
    # Calculate intermediate values as in the proof generation
    mu = alpha * gamma

    # Create delta vector
    delta_vec = vector(Fq, n)
    for i in range(r):
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

# =============================================================================
# GENERATOR SETUP (WILL BE INITIALIZED BY SPECIFIC IMPLEMENTATIONS)
# =============================================================================

# These will be set by the specific eVRF implementations
G_vec = None
H_vec = None
G = None
H = None
GT1 = None
GT2 = None
GS = None