import hashlib

# base fild of BLS12-381
p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
# order of BLS12-381
r = 52435875175126190479447740508185965837690552500527637822603658699938581184513
# order bandersnatch
q = 15172417585395309745210573063711216967055694857434315578142854216712503379

assert p in Primes()
assert r in Primes()
assert q in Primes()

Fp = GF(p)
Fr = GF(r)
Fq = GF(q)

# BLS12-381 G1
E = EllipticCurve(Fp, [0, 4])
# Bandersnatch
F = EllipticCurve_from_j(Fr(8000))
a = 52435875175126190479447740508185965837690552500527637822603658699938430656513
b = 52435875175126190479447740508185965837690552500527637822603658699309173440513

# G1 generator
G = E([0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb,
0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1])
assert (r*G).is_zero()

#private key 
#k =  Fq.random_element()
  
# this is the naive hash and pray
# TODO implement a proper hash_to_curve
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

def keygen():

    #Private key
    k = Fr.random_element()

    #verification key
    Q = k * G

    #Schnorr proof to generate pi_Q
    pi_Q = proof_R_dlog(k)

    vk = (Q, pi_Q)

    return k, vk

#Schnorr proof of knowledge of discrete log
def proof_R_dlog(k):

    k_fr = Fr(k)
    P = k_fr*G

    r_nonce = Fr.random_element() #Random nonce
    R = r_nonce*G 

    challenge_data = point_to_bytes(P) + point_to_bytes(G) + point_to_bytes(R)

    challenge = Fr(int.from_bytes(hashlib.sha256(challenge_data).digest(), byteorder='big'))

    #Challeng response
    s = r_nonce + challenge * k_fr



    return (R, s)


#Verify Discrete log proof that Q = k*G for some k
def verify_proof_R_dlog(P, proof):
     
    R, s = proof

    #Compute the challenge as in the proof

    challenge_data = point_to_bytes(P) + point_to_bytes(G) + point_to_bytes(R)

    challenge = Fr(int.from_bytes(hashlib.sha256(challenge_data).digest(), byteorder='big'))

    left = Integer(s)*G
    right = R + challenge * P

    return left == right


# Generate a complete witness vector for the R1CS constraint system (equation 12)
def generate_witness_vector(k, H_x, G_S):
    ell = 256
    
    # Size of the witness vector
    m = 8*ell + 6
    
    # Initialize the witness vector
    z = vector(Fr, m)
    
    # First element is the constant 1
    z[0] = 1
    
    # Second element is the private key k
    z[1] = k
    
    # Third element is the x-coordinate of the final point P = k*H_x
    final_point = k * H_x
    z[2] = final_point[0]  
    
    # Define the constants c_i as per equation (9)
    c = [i+2 for i in range(ell - 1)]
    c.append(-sum(c))  # Make the sum zero
    
    # Convert k to its binary representation (LSB first)
    k_bits = Integer(k).digits(base=2, padto=ell)
    
    # Compute all points P_i according to equation (10)
    P = []
    
    P.append(k_bits[0] * H_x + c[0] * G_S)
    
    for i in range(1, ell):

        delta_i = (2^i * k_bits[i]) * H_x + c[i] * G_S
        P.append(P[i-1] + delta_i)


    # append
    z[3] = k_bits[0]

    z[4] = P[0][0]

    z[5] = P[0][1]
    
    print("length of p:", len(P))

    for i in range(1, ell):
        base_idx = 8*i - 2
        
        z[base_idx] = k_bits[i]
        
        z[base_idx + 1] = P[i][0]
        
        z[base_idx + 2] = P[i][0]^2
        
        z[base_idx + 3] = P[i][0]^3
        
        z[base_idx + 4] = P[i][1]
        
        z[base_idx + 5] = P[i][1]^2
        
        # For each bit, we need to compute t_1 and t_2
        
        prev_point = P[i-1]
        
        delta_i = (2^i * k_bits[i]) * H_x + c[i] * G_S
        delta_prime_i = c[i] * G_S
        
        x_prev = prev_point[0]
        y_prev = prev_point[1]
        x_curr = P[i][0]
        y_curr = P[i][1]
        
        delta_x = delta_i[0] - delta_prime_i[0]
        delta_y = delta_i[1] - delta_prime_i[1]
        x_prime = delta_prime_i[0]
        y_prime = delta_prime_i[1]
        
        # Compute t_1 = (δ_x*k_i + x′ - x_P_i)*(y_P_i + y_P_{i-1})
        t_1 = (delta_x*k_bits[i] + x_prime - x_curr) * (y_curr + y_prev)
        
        
        # Compute t_2 = (x_P_{i-1} - x_P_i)*(δ_y*k_i + y′ + y_P_i)
        t_2 = (x_prev - x_curr) * (delta_y*k_bits[i] + y_prime + y_curr)
        

        # Add t_1 and t_2 to the witness vector
        z[base_idx + 6] = t_1
        z[base_idx + 7] = t_2
        
    
    return z

# Construct the R1CS matrices A, B, C for the DDH relation
def construct_R1CS_matrices(H_x, G_S):
    # Use 256 bits for the key length
    ell = 256
    
    # Dimensions of the matrices with optimization
    n = 5*ell + 2  
    m = 8*ell + 6  
    
    # Initialize matrices with zeros
    A = matrix(Fr, n, m)
    B = matrix(Fr, n, m)
    C = matrix(Fr, n, m)
    
    # Define constants c_i as per equation (9) in the paper
    c = [i+2 for i in range(ell-1)]
    c.append(-sum(c))  # Make the sum zero
    
    # For each bit position, compute the delta values
    delta_x = []
    delta_y = []
    x_prime = []
    y_prime = []
    
    for i in range(ell):
        # Compute the point Delta_i = (2^i) * H_x + c_i * G_S
        Delta_full = (2^i) * H_x + c[i] * G_S
        Delta_x = Delta_full[0]
        Delta_y = Delta_full[1]
        
        # The point Delta'_i = c_i * G_S
        Delta_prime_full = c[i] * G_S
        Delta_prime_x = Delta_prime_full[0]
        Delta_prime_y = Delta_prime_full[1]
        
        # Compute deltas
        delta_x.append(Delta_x - Delta_prime_x)
        delta_y.append(Delta_y - Delta_prime_y)
        x_prime.append(Delta_prime_x)
        y_prime.append(Delta_prime_y)
    
    # Current constraint index
    constraint_idx = 0
    
    # First constraint: key = sum(2^i * k_i)
    A[constraint_idx, 3] = 2^0  # 2^0 * k_0
    for i in range(1, ell):
        A[constraint_idx, 8*i-2] = 2^i  # 2^i * k_i
    
    B[constraint_idx, 0] = 1  # 1
    C[constraint_idx, 1] = 1  # k
    constraint_idx += 1
    
    # Add constraint to check that the final point P_ℓ has the correct x-coordinate
    A[constraint_idx, 8*(ell-1)-1] = 1  # x_P_ℓ (the x-coordinate of the final point)
    B[constraint_idx, 0] = 1  # 1
    C[constraint_idx, 2] = 1  # x_P (the output)
    constraint_idx += 1
    
    # For each bit of the key, add constraints
    for i in range(ell):
        if i == 0:
            k_idx = 3  # k_0 is at z[3]
            x_idx = 4  # x_P_0 is at z[4]
            y_idx = 5  # y_P_0 is at z[5]
            # Note: For i=0, we don't have x^2, x^3, y^2, t1, t2 in the witness
        else:
            k_idx = 8*i-2  # k_i
            x_idx = 8*i-1  # x_P_i
            x2_idx = 8*i   # x_P_i^2
            x3_idx = 8*i+1  # x_P_i^3
            y_idx = 8*i+2  # y_P_i
            y2_idx = 8*i+3  # y_P_i^2
            t1_idx = 8*i+4  # t1
            t2_idx = 8*i+5  # t2
        
        # Previous point indices
        if i > 0:
            x_prev_idx = 8*(i-1)-1 if i > 1 else 4  # x_P_{i-1}
            y_prev_idx = 8*(i-1)+2 if i > 1 else 5  # y_P_{i-1}
        else:
            x_prev_idx = m-2  # For i=0, use dummy indices
            y_prev_idx = m-1
        
        # 1. Constraint for k_i * (1 - k_i) = 0 (enforcing k_i is binary)
        A[constraint_idx, k_idx] = 1  # k_i
        B[constraint_idx, 0] = 1  # 1 (constant)
        B[constraint_idx, k_idx] = -1  # -k_i
        C[constraint_idx, 0] = 0  # Result should be 0
        constraint_idx += 1
        
        # Skip remaining constraints for i=0 since we don't have those values in witness
        if i == 0:
            continue
        
        # Combined constraint for x squaring and cubing
        alpha = Fr.random_element()
        
        A[constraint_idx, x_idx] = 1  # x_P_i
        B[constraint_idx, x_idx] = 1  # x_P_i
        C[constraint_idx, x2_idx] = 1  # x_P_i^2
        
        # Add the second constraint with weight alpha
        A[constraint_idx, x_idx] += alpha  # alpha * x_P_i
        B[constraint_idx, x2_idx] += alpha  # alpha * x_P_i^2
        C[constraint_idx, x3_idx] = alpha  # alpha * x_P_i^3
        
        constraint_idx += 1
        
        # Constraint for y_P_i * y_P_i = y_P_i^2
        A[constraint_idx, y_idx] = 1  # y_P_i
        B[constraint_idx, y_idx] = 1  # y_P_i
        C[constraint_idx, y2_idx] = 1  # y_P_i^2
        constraint_idx += 1
        
        # Combined constraint for elliptic curve equation and t_1 = t_2
        beta = Fr.random_element()
        
        # Elliptic curve equation component
        A[constraint_idx, x3_idx] = 1  # x_P_i^3
        A[constraint_idx, x_idx] = a  # a*x_P_i
        A[constraint_idx, 0] = b  # b (constant)
        B[constraint_idx, 0] = 1  # 1 (for the constant term)
        C[constraint_idx, y2_idx] = 1  # y_P_i^2
        
        # t_1 = t_2 component with weight beta
        A[constraint_idx, t1_idx] += beta  # beta * t_1
        B[constraint_idx, 0] += 0  # No change needed
        C[constraint_idx, t2_idx] += beta  # beta * t_2
        
        constraint_idx += 1
        
        # Combined constraint for both slope computations
        gamma = Fr.random_element()
        
        # First slope constraint component
        A[constraint_idx, k_idx] = delta_x[i]  # δ_x*k_i
        A[constraint_idx, 0] = x_prime[i]  # x′
        A[constraint_idx, x_idx] = -1  # -x_P_i
        B[constraint_idx, y_idx] = 1  # y_P_i
        B[constraint_idx, y_prev_idx] = 1  # y_P_{i-1}
        C[constraint_idx, t1_idx] = 1  # t_1
        
        # Second slope constraint component with weight gamma
        A[constraint_idx, x_prev_idx] += gamma  # gamma * x_P_{i-1}
        A[constraint_idx, x_idx] += -gamma  # gamma * (-x_P_i)
        B[constraint_idx, k_idx] += gamma * delta_y[i]  # gamma * δ_y*k_i
        B[constraint_idx, 0] += gamma * y_prime[i]  # gamma * y′
        B[constraint_idx, y_idx] += gamma  # gamma * y_P_i
        C[constraint_idx, t2_idx] = gamma  # gamma * t_2
        
        constraint_idx += 1
    
    assert constraint_idx <= n, f"Expected at most {n} constraints, but created {constraint_idx}"
    
    return A, B, C





def proof_BP(A, B, C, T, z):
    
    # Step 1: Extract dimensions
    n, m = A.nrows(), A.ncols()
    
    # Step 2: Compute Az, Bz, Cz to verify the constraints
    Az = A * z
    Bz = B * z
    Cz = C * z

    # H is generator in source group
    H = F.random_point()
    
    
    # Verify that z satisfies the R1CS relationship (Az) ◦ (Bz) = (Cz)
    for i in range(n):
        if Az[i] * Bz[i] != Cz[i]:
           raise ValueError(f"Witness does not satisfy R1CS constraint {i}")
    
    # Initialize the Bulletproofs protocol
    
    # Generate random blinding factors for Pedersen commitments
    r_a = Fr.random_element()
    r_b = Fr.random_element()
    r_c = Fr.random_element()
    
    # Compute Pedersen commitments to Az, Bz, Cz
    com_a = sum(Az[i] * G for i in range(n)) + r_a * H
    com_b = sum(Bz[i] * G for i in range(n)) + r_b * H
    com_c = sum(Cz[i] * G for i in range(n)) + r_c * H
    
    hash_input = str(T) + str(com_a) + str(com_b)
    y_seed = int.from_bytes(hashlib.sha256(hash_input.encode()).digest(), byteorder='big')
    y = Fr(y_seed % Fr.order())
    
    # Compute the combined witness vectors
    a_combined = vector(Fr, [Az[i] * y^i for i in range(n)])
    b_combined = vector(Fr, [Bz[i] * y^i for i in range(n)])
    c_combined = vector(Fr, [Cz[i] * y^i for i in range(n)])
    
    # Generate random blinding factors for the inner product argument
    alpha = Fr.random_element()
    beta = Fr.random_element()
    
    # Generate vectors for inner product
    a_vec = vector(Fr, [a_combined[i] for i in range(n)])
    b_vec = vector(Fr, [b_combined[i] for i in range(n)])
    
    # Perform the logarithmic-sized inner product proof
    L_values = []
    R_values = []
    challenges = []
    
    # Number of rounds needed for the inner product argument
    num_rounds = (n-1).nbits()
    
    # Current vectors in the protocol
    a_current = a_vec
    b_current = b_vec
    
    for i in range(num_rounds):
        # Split vectors in half
        n_half = len(a_current) // 2
        a_left = a_current[:n_half]
        a_right = a_current[n_half:]
        b_left = b_current[:n_half]
        b_right = b_current[n_half:]
        
        # Compute inner products
        c_left = sum(a_right[j] * b_left[j] for j in range(len(a_right)))
        c_right = sum(a_left[j] * b_right[j] for j in range(len(a_left)))
        
        # Compute L and R values (Pedersen commitments)
        r_L = Fr.random_element()
        r_R = Fr.random_element()
        
        L = c_left * G + r_L * H
        R = c_right * G + r_R * H
        
        L_values.append(L)
        R_values.append(R)
        
        # Generate challenge for this round
        x = Fr.random_element()
        challenges.append(x)
        
        # Update the vectors for the next round
        a_new = vector(Fr, [a_left[j] * x + a_right[j] / x for j in range(n_half)])
        b_new = vector(Fr, [b_left[j] / x + b_right[j] * x for j in range(n_half)])
        
        a_current = a_new
        b_current = b_new
    
    
    final_a = a_current[0]
    final_b = b_current[0]
    
    # Assemble the Bulletproof
    proof = {
        "L_values": L_values,
        "R_values": R_values,
        "final_a": final_a,
        "final_b": final_b
    }
    
    return proof

def verify_BP(A, B, C, T, proof):
    #TODO: Implement the verification of the Bulletproof
    # Extract the proof components
   

    return True

def proof_relation_H(Q, x, Y, k):
    
    H_x = hash_to_curve_bandersnatch(x)
    
    A, B, C = construct_R1CS_matrices(H_x, F.random_point())
    
    T = Q + Y
    
    z = generate_witness_vector(k, H_x, F.random_point())
    
    # Generate the Bulletproof
    pi_0 = proof_BP(A, B, C, T, z)
    
    # Generate the discrete log proofs
    pi_Q = proof_R_dlog(k)
    pi_Y = proof_R_dlog(Y)
    
    return (pi_0, pi_Q, pi_Y)

def verify_proof_relation_H(Q, x, Y, pi):

    pi_0, pi_Q, pi_Y = pi

    H_x = hash_to_curve_bandersnatch(x)


    #Verify that Q has a valid proof
    if not verify_proof_R_dlog(Q, pi_Q):
        return False
    
    if not verify_proof_R_dlog(Y, pi_Y):
        return False

    # Construct the R1CS matrices
    A, B, C = construct_R1CS_matrices(H_x, G)

    T = Q + Y 
    
    if not verify_BP(A, B, C, T, pi_0):
        return False

    return True

def eval(k,x):
    H_x = hash_to_curve_bandersnatch(x)

    P = k * H_x
    
    # Extract the x-coordinate of P
    x_P = P[0]
    
   
    y = x_P
    
    # Compute Y = y·G
    Y = y * G
    
    # Compute Q = k·G
    Q = k * G
    
    pi = proof_relation_H(Q, x, Y, k)
    
    #Generate proof of discrete log of y
    #pi = proof_R_dlog(y)

    return y, Y, pi

def verify(vk, x, Y, pi):

    Q, pi_Q = vk

    return verify_proof_relation_H(Q, x, Y, pi)

""" def test():
    #Generate key
    k, vk = keygen()

    #print("Generated key pair:")
    #print(f"Secret key k: {k}")
    #print(f"Verification key Q: {vk[0]}")

    # Evaluate the eVRF on some input
    x = Fr.random_element()
    y, Y, pi = eval(k, x)
    #print("\nEvaluated eVRF on input x =", x)
    #print(f"Output y: {y}")
    #print(f"Output Y: {Y}")

    # Verify the eVRF output
    result = verify(vk, x, Y, pi)
    print("\nVerification result:", result)

test()
"""


# Test cases for the EVRF implementation
def test_keygen():
    k, vk = keygen()
    Q, pi_Q = vk
    
    # Test that k is in Fr
    assert k in Fr
    
    # Test that Q is on curve E
    assert Q in E
    
    # Test that Q = k*G
    assert Q == k*G
    
    # Test that the proof verifies
    assert verify_proof_R_dlog(Q, pi_Q)

def test_hash_to_curve():
    x = Fr.random_element()
    P = hash_to_curve_bandersnatch(x)
    
    # Test that output point is on curve F
    assert P in F

def test_proof_dlog():
    # Test with known discrete log
    k = Fr.random_element()
    P = k * G
    
    proof = proof_R_dlog(k)
    
    # Verify the proof
    assert verify_proof_R_dlog(P, proof)
    
    # Test with wrong point
    wrong_P = (k+_sage_const_1 ) * G
    assert not verify_proof_R_dlog(wrong_P, proof)

def test_witness_vector():
    k = Fr.random_element()
    x = Fr.random_element()
    H_x = hash_to_curve_bandersnatch(x)
    G_S = F.random_point()
    
    z = generate_witness_vector(k, H_x, G_S)

    
    
    # Test dimensions
    assert len(z) == 2054
    
    # Test first elements
    assert z[0] == 1 
    assert z[1] == k
    assert z[2] == (k * H_x)[0]

def test_R1CS_matrices():
    x = Fr.random_element()
    H_x = hash_to_curve_bandersnatch(x)
    G_S = F.random_point()
    
    A, B, C = construct_R1CS_matrices(H_x, G_S)
    
    # Test dimensions
    n = 5 *256  + 2 
    m = 8 *256  + 6 
    
    assert A.dimensions() == (n, m)
    assert B.dimensions() == (n, m)
    assert C.dimensions() == (n, m)

def test_eval_verify():
    # Generate keys
    k, vk = keygen()
    
    # Test evaluation
    x = Fr.random_element()
    y, Y, pi = eval(k, x)
    
    # Test verification
    assert verify(vk, x, Y, pi)
    
    # Test with wrong input
    wrong_x = x + 1 
    assert not verify(vk, wrong_x, Y, pi)
    
    # Test with wrong output point
    wrong_Y = Y + G
    assert not verify(vk, x, wrong_Y, pi)

def run_all_tests():
    print("Running tests...")
    
    try:
        test_keygen()
        print("✓ keygen test passed")
        
        test_hash_to_curve()
        print("✓ hash_to_curve test passed")
        
        test_proof_dlog()
        print("✓ discrete log proof test passed")
        
        test_witness_vector()
        print("✓ witness vector test passed")
        
        test_R1CS_matrices()
        print("✓ R1CS matrices test passed")
        
        test_eval_verify()
        print("✓ eval/verify test passed")
        
        print("\nAll tests passed successfully!")
        
    except AssertionError as e:
        print(f"Test failed: {str(e)}")
    except Exception as e:
        print(f"Unexpected error: {str(e)}")

if __name__ == "__main__":
    run_all_tests()




