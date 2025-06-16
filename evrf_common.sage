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
    """Generate modified inner product proof (recursive) following the paper"""
    n = len(a)
    if n == 1:
        # Base case: vectors of length 1
        r = Fq.random_element()
        s = Fq.random_element()
        delta = Fq.random_element()
        eta = Fq.random_element()

       
        # Compute commitments A and B
        A = r * G[0] + s * H[0] + (r * b[0] + s * a[0]) * G_base + delta * H_base
        B = (r * s) * G_base + eta * H_base
        
        # Generate Fiat-Shamir challenge
        challenge_data = point_to_bytes(A) + point_to_bytes(B)
        e = Fq(int.from_bytes(hashlib.sha256(challenge_data).digest(), byteorder='big'))
        
        P_computed = inner_product(a, G) + inner_product(b, H) + inner_product(a, b) * G_base + alpha * H_base


        # Compute responses
        r_prime = r + a[0] * e
        s_prime = s + b[0] * e
        delta_prime = eta + delta * e + alpha * e^2

        
        

        
        
        return {
            "A": A,
            "B": B,
            "e": e,
            "P_computed": P_computed,
            "r_prime": r_prime,
            "s_prime": s_prime,
            "delta_prime": delta_prime
        }
    
    # Recursive case
    n_half = n // 2
    
    # Split vectors
    a_1 = vector(Fq, a[:n_half]) 
    a_2 = vector(Fq, a[n_half:])  
    b_1 = vector(Fq, b[:n_half])  
    b_2 = vector(Fq, b[n_half:])  
    G_1 = G[:n_half]               
    G_2 = G[n_half:]               
    H_1 = H[:n_half]               
    H_2 = H[n_half:]               
    
    # Compute inner products for cross terms
    c_L = inner_product(a_1, b_2)  
    c_R = inner_product(a_2, b_1)  
    
    # Sample random blinding factors
    d_L = Fq.random_element()
    d_R = Fq.random_element()
    
    # Compute commitments L and R according to the paper
    # L = <a_1, G_2> + <b_2, H_1> + c_L * G + d_L * H
    L = inner_product(a_1, G_2) + inner_product(b_2, H_1) + c_L * G_base + d_L * H_base
    
    # R = <a_2, G_1> + <b_1, H_2> + c_R * G + d_R * H
    R = inner_product(a_2, G_1) + inner_product(b_1, H_2) + c_R * G_base + d_R * H_base
    
    # Generate Fiat-Shamir challenge
    challenge_data = point_to_bytes(L) + point_to_bytes(R)
    e = Fq(int.from_bytes(hashlib.sha256(challenge_data).digest(), byteorder='big'))
    e_inv = e^(-1)
    
    
    # G' = e^(-1) * G_1 + e * G_2
    # H' = e * H_1 + e^(-1) * H_2
    G_prime = []
    H_prime = []
    for i in range(n_half):
        G_prime.append(e_inv * G_1[i] + e * G_2[i])
        H_prime.append(e * H_1[i] + e_inv * H_2[i])
    
    # Compute new vectors for the next level
    # a' = e * a_1 + e^(-1) * a_2
    # b' = e^(-1) * b_1 + e * b_2
    a_prime = []
    b_prime = []
    for i in range(n_half):
        a_prime.append(e * a_1[i] + e_inv * a_2[i])
        b_prime.append(e_inv * b_1[i] + e * b_2[i])
    
    
    P_prime = e^2 * L + P + e^(-2) * R


    alpha_prime = e^2 * d_L + alpha + e^(-2) * d_R
    
    # Recursively generate the inner product proof
    inner_proof = generate_modified_inner_product_proof(
        G_prime, H_prime, G_base, H_base, P_prime, 
        a_prime, b_prime, alpha_prime
    )
    
    
    return {
        "L": L,
        "R": R,
        "e": e,
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
        
        A = proof["A"]
        B = proof["B"]
        e = proof["e"]
        P_computed = proof["P_computed"]
        r_prime = proof["r_prime"]
        s_prime = proof["s_prime"]
        delta_prime = proof["delta_prime"]
        
        # Verify the challenge matches what we'd compute
        challenge_data = point_to_bytes(A) + point_to_bytes(B)
        e_computed = Fq(int.from_bytes(hashlib.sha256(challenge_data).digest(), byteorder='big'))

        if e != e_computed:
            return False
        
        left_side = e^2 * P + e * A + B

        left_side_with_P = e^2 * P_computed + e * A + B
        
        # Build right side term by term
        term1 = (r_prime * e )* G[0]         
        term2 = (s_prime * e )* H[0]         
        term3 = (r_prime * s_prime) * G_base 
        term4 = delta_prime * H_base       
        
        right_side = term1 + term2 + term3 + term4

        print("Left side:", left_side)
        print("Right side:", right_side)
        print("Left side with P_computed:", left_side_with_P)
        
        
        return left_side == right_side
    
    # Recursive case
    
    L = proof["L"]
    R = proof["R"]
    e = proof["e"]
    inner_proof = proof["inner_proof"]
    
    
    challenge_data = point_to_bytes(L) + point_to_bytes(R)
    e_computed = Fq(int.from_bytes(hashlib.sha256(challenge_data).digest(), byteorder='big'))
    if e != e_computed:
        return False
    
    # Calculate the new values for the recursive verification
    n_half = n // 2
    G_1 = G[:n_half]
    G_2 = G[n_half:]
    H_1 = H[:n_half]
    H_2 = H[n_half:]
    
    e_inv = e^(-1)

    # Compute new generators for the next level
    # G' = e^(-1) * G_1 + e * G_2
    # H' = e * H_1 + e^(-1) * H_2
    G_prime = []
    H_prime = []
    for i in range(n_half):
        G_prime.append(e_inv * G_1[i] + e * G_2[i])
        H_prime.append(e * H_1[i] + e_inv * H_2[i])
    
    # Compute the new point P' = e^2 * L + P + e^(-2) * R
    P_prime = e^2 * L + P + e^(-2) * R
    
    
    # Recursively verify the inner product proof
    return verify_modified_inner_product_proof(G_prime, H_prime, G_base, H_base, P_prime, inner_proof)


def generate_bulletproof(A, B, C, T, z):
    """Generate bulletproof for R1CS statement following the paper's R1CS* protocol"""
    m_b, n_b = A.nrows(), A.ncols()

    # In this implementation, z' is a zero vector and eta is 0 in Fq
    z_prime = vector(Fq, n_b, [0] * n_b)
    eta = Fq(0)

    # r is the size of x in the witness vector z = (x||y)
    r = 3
    x = vector(Fq, [z[i] for i in range(r)])
    y = vector(Fq, [z[i] for i in range(r, n_b)])
    
    # Sample random r for S
    r_random = Fq.random_element()
    
    # Compute Az and Bz
    Az = A * z
    Bz = B * z
    
    # Step 1: Prover computes S
    # S = ⟨((x'||y) || Az), G⟩ + ⟨(0^n || Bz), H⟩ + r·H
    
    x_prime = vector(Fq, [0] * r) 
    x_prime_y_concat = vector(Fq, list(x_prime) + list(y))  
    
    # First part: ⟨((x'||y) || Az), G⟩
    G_term_vec = vector(Fq, list(x_prime_y_concat) + list(Az))  
    S_G_term = inner_product(G_term_vec, G_vec) 
    
    # Second part: ⟨(0^n || Bz), H⟩  
    zero_n = vector(Fq, [0] * n_b)
    H_term_vec = vector(Fq, list(zero_n) + list(Bz))  
    S_H_term = inner_product(H_term_vec, H_vec)
    
    # Complete S computation
    S = S_G_term + S_H_term + r_random * H
    
    # Step 2: Use Fiat-Shamir to get verifier's challenges
    transcript_data = point_to_bytes(S)
    alpha = compute_challenge([transcript_data, b"alpha"])
    beta = compute_challenge([transcript_data, b"beta"])
    gamma = compute_challenge([transcript_data, b"gamma"])
    delta = compute_challenge([transcript_data, b"delta"])
    
    # Step 3: Compute intermediate values as per the paper
    mu = alpha * gamma
    
    
    delta_vec = vector(Fq, n_b)
    for i in range(r):
        delta_vec[i] = delta
    for i in range(r, n_b):
        delta_vec[i] = 1
    
    
    # Construct δ^(-1) vector
    delta_inv = delta^(-1)
    delta_inv_vec = vector(Fq, n_b)
    for i in range(r):
        delta_inv_vec[i] = delta_inv
    for i in range(r, n_b):
        delta_inv_vec[i] = 1
    
    
    # Step 4: Construct G' vector with modified generators
    # G' = (G_1, ..., G_n, γ^(-1)·G_{n+1}, ..., γ^(-m)·G_{n+m})
    gamma_inv = gamma^(-1)
    G_prime = []
    
    for i in range(n_b):
        G_prime.append(G_vec[i])
    
    for i in range(m_b):
        G_prime.append((gamma_inv^(i+1)) * G_vec[n_b + i])

    
    # Step 5: Compute vector c = μ^m A + β^m B - γ^m C
    c = vector(Fq, n_b)
    
    mu_vec = vector(Fq, [mu^(i+1) for i in range(m_b)])
    beta_m_vec = vector(Fq, [beta^(i+1) for i in range(m_b)])
    gamma_m_vec = vector(Fq, [gamma^(i+1) for i in range(m_b)])

    c = mu_vec * A + beta_m_vec * B - gamma_m_vec * C
    
    assert len(c) == n_b, "Vector c must have the same length as the number of columns in A, B, C"
    # Step 6: Compute inner product value ω
    # ω = ⟨α^m, β^m⟩ + δ² · ⟨α^n, c ◦ δ⟩

    alpha_m_vec = vector(Fq, [alpha^(i+1) for i in range(m_b)])
    omega = 0
    
    omega += inner_product(alpha_m_vec, beta_m_vec)  # First term: ⟨α^m, β^m⟩
    
    alpha_n_vec = vector(Fq, n_b, [alpha^(i+1) for i in range(n_b)])
    
    c_hadamard_delta = hadamard_product(c, delta_vec)
    omega += (delta^2) * inner_product(alpha_n_vec, c_hadamard_delta)
    
    
    P = delta_inv * T + S
    
    # Add ⟨(δ² · α^n || -β^m), G'⟩
    delta_sq_alpha_n = (delta^2) * alpha_n_vec
    neg_beta_m = - beta_m_vec
    g_prime_term_vec = vector(Fq, list(delta_sq_alpha_n) + list(neg_beta_m))
    P += inner_product(g_prime_term_vec, G_prime)
    
    # Add ⟨(c ◦ δ || -α^m), H⟩
    neg_alpha_m = - alpha_m_vec
    h_term_vec = vector(Fq, list(c_hadamard_delta) + list(neg_alpha_m))
    P += inner_product(h_term_vec, H_vec)

    
    # For u: ((x'||y) + δ^(-1) · (x||y') + δ² · α^n || (Az + δ^(-1) · Az') ◦ γ^m - β^m)
    y_prime = vector(Fq, [0] * (n_b - r)) 
    x_y_prime_concat = vector(Fq, list(x) + list(y_prime))
    
    # First part of u
    u_first_part = x_prime_y_concat + delta_inv * x_y_prime_concat + (delta^2) * alpha_n_vec
    
    # Second part of u
    Az_prime = A * z_prime
    Az_combined = Az + delta_inv * Az_prime
    gamma_m_vec = vector(Fq, [gamma^(i+1) for i in range(m_b)])
    Az_gamma_hadamard = hadamard_product(Az_combined, gamma_m_vec)
    u_second_part = Az_gamma_hadamard - beta_m_vec
    
    u = vector(Fq, list(u_first_part) + list(u_second_part))
    
    # For v: (c ◦ δ || Bz - α^m + δ^(-1) · Bz')
    Bz_prime = B * z_prime
    v_second_part = Bz - alpha_m_vec + delta_inv * Bz_prime
    
    v = vector(Fq, list(c_hadamard_delta) + list(v_second_part))
    
    # Step 8: Compute η'
    eta_prime = r_random + delta_inv * eta  

    P_computed = inner_product(u, G_prime) + inner_product(v, H_vec) + eta_prime * H

    print("P_computed:", P_computed)
    print("P:", P)

   
    # Step 11: Generate the inner product proof
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
    """Verify bulletproof for R1CS statement following the paper's R1CS* protocol"""
    # Extract proof components
    S = proof["S"]
    challenges = proof["challenges"]
    alpha = challenges["alpha"]
    beta = challenges["beta"]
    gamma = challenges["gamma"]
    delta = challenges["delta"]
    inner_product_proof = proof["inner_product_proof"]

    m_b, n_b = A.nrows(), A.ncols()
    r = 3  # Size of x in the witness vector z = (x||y)

    mu = alpha * gamma

    delta_vec = vector(Fq, n_b)
    for i in range(r):
        delta_vec[i] = delta
    for i in range(r, n_b):
        delta_vec[i] = 1
    
    
    # Construct δ^(-1) vector
    delta_inv = delta^(-1)
    delta_inv_vec = vector(Fq, n_b)
    for i in range(r):
        delta_inv_vec[i] = delta_inv
    for i in range(r, n_b):
        delta_inv_vec[i] = 1
    
    
    # Step 4: Construct G' vector with modified generators
    # G' = (G_1, ..., G_n, γ^(-1)·G_{n+1}, ..., γ^(-m)·G_{n+m})
    gamma_inv = gamma^(-1)
    G_prime = []
    
    for i in range(n_b):
        G_prime.append(G_vec[i])
    
    for i in range(m_b):
        G_prime.append((gamma_inv^(i+1)) * G_vec[n_b + i])

    
    # Step 5: Compute vector c = μ^m A + β^m B - γ^m C
    c = vector(Fq, n_b)
    
    mu_vec = vector(Fq, [mu^(i+1) for i in range(m_b)])
    beta_m_vec = vector(Fq, [beta^(i+1) for i in range(m_b)])
    gamma_m_vec = vector(Fq, [gamma^(i+1) for i in range(m_b)])

    c = mu_vec * A + beta_m_vec * B - gamma_m_vec * C
    
    assert len(c) == n_b, "Vector c must have the same length as the number of columns in A, B, C"
    # Step 6: Compute inner product value ω
    # ω = ⟨α^m, β^m⟩ + δ² · ⟨α^n, c ◦ δ⟩

    alpha_m_vec = vector(Fq, [alpha^(i+1) for i in range(m_b)])
    omega = 0
    
    omega += inner_product(alpha_m_vec, beta_m_vec)  # First term: ⟨α^m, β^m⟩
    
    alpha_n_vec = vector(Fq, n_b, [alpha^(i+1) for i in range(n_b)])
    
    c_hadamard_delta = hadamard_product(c, delta_vec)
    omega += (delta^2) * inner_product(alpha_n_vec, c_hadamard_delta)
    
    
    P = delta_inv * T + S
    
    # Add ⟨(δ² · α^n || -β^m), G'⟩
    delta_sq_alpha_n = (delta^2) * alpha_n_vec
    neg_beta_m = - beta_m_vec
    g_prime_term_vec = vector(Fq, list(delta_sq_alpha_n) + list(neg_beta_m))
    P += inner_product(g_prime_term_vec, G_prime)
    
    # Add ⟨(c ◦ δ || -α^m), H⟩
    neg_alpha_m = - alpha_m_vec
    h_term_vec = vector(Fq, list(c_hadamard_delta) + list(neg_alpha_m))
    P += inner_product(h_term_vec, H_vec)

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