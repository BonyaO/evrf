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

    print("Values in proof")
    print("left")
    print(s*G)
    print(R + challenge * P)

    return (R, s)


#Verify Discrete log proof that Q = k*G for some k
def verify_proof_R_dlog(P, proof):
     
    R, s = proof

    #Compute the challenge as in the proof

    challenge_data = point_to_bytes(P) + point_to_bytes(G) + point_to_bytes(R)

    challenge = Fr(int.from_bytes(hashlib.sha256(challenge_data).digest(), byteorder='big'))

    left = Integer(s)*G
    right = R + challenge * P

    print("Values in Verification")
    print("left")
    print(left)
    print("right")
    print(right)

    return left == right


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
    
    # Generate witness vector
    # z = generate_witness_vector(k, H_x, F.random_point())
    
    # Generate R1CS matrices
    # A, B, C = construct_R1CS_matrices(H_x, F.random_point())

    #TODO: Generate bulletproof argument
    # pi = generate_bulletproof_argument(A, B, C, z, G, G) 


    return y, Y

def verify(vk, x, Y):

    Q, pi_Q = vk

    #Verify that Q has a valid proof
    if not verify_proof_R_dlog(Q, pi_Q):
        return False
    

    #TODO: Implemtne the verify_proof_R_H function
    #result = verify_proof_R_H(Q, x, Y, pi)

    return True

def test():
    #Generate key
    k, vk = keygen()

    #print("Generated key pair:")
    #print(f"Secret key k: {k}")
    #print(f"Verification key Q: {vk[0]}")

    # Evaluate the eVRF on some input
    x = Fr.random_element()
    y, Y = eval(k, x)
    #print("\nEvaluated eVRF on input x =", x)
    #print(f"Output y: {y}")
    #print(f"Output Y: {Y}")

    # Verify the eVRF output
    result = verify(vk, x, Y)
    print("\nVerification result:", result)

test()


# Construct the R1CS matrices A, B, C for the DDH relation
def construct_R1CS_matrices(H_x, G_S):
    
    # Use 256 bits for the key length
    ell = 256
    
    # Dimensions of the matrices
    n = 5*ell + 2  # number of constraints
    m = 8*ell + 6  # size of the witness vector
    
    
    A = matrix(Fq, n, m)
    B = matrix(Fq, n, m)
    C = matrix(Fq, n, m)
    
    # Get curve parameters
    # Bandersnatch curve in Weierstrass form: y^2 = x^3 + ax + b
    a = 52435875175126190479447740508185965837690552500527637822603658699938430656513
    b = 52435875175126190479447740508185965837690552500527637822603658699309173440513
    
    # Define the constants c_i as per equation (9)
    c = [i+2 for i in range(ell)]
    c.append(-sum(c))  # Make the sum zero
    
    # Define delta_x and delta_y for each bit
    delta_x = []
    delta_y = []
    x_prime = []
    y_prime = []
    
    # For each bit, compute the specific values
    for i in range(ell):
        # Compute the point Delta_i = (2^i) * H_x + c_i * G_S
        Delta_full = (2^i) * H_x + c[i] * G_S
        Delta_x = Delta_full[0]
        Delta_y = Delta_full[1]
        
        # The point Delta'_i = c_i * G_S
        Delta_prime_full = c[i] * G_S
        Delta_prime_x = Delta_prime_full[0]
        Delta_prime_y = Delta_prime_full[1]
        
        # Compute deltas as in equation (13)
        delta_x.append(Delta_x - Delta_prime_x)
        delta_y.append(Delta_y - Delta_prime_y)
        x_prime.append(Delta_prime_x)
        y_prime.append(Delta_prime_y)
    
    # Current constraint index
    constraint_idx = 0
    
    # For each bit of the key, add the constraints
    for i in range(ell):
        # Indices for the witness elements for this bit
        k_idx = 2 + 8*i         # k_i
        x_idx = 3 + 8*i         # x_P_i
        x2_idx = 4 + 8*i        # x_P_i^2
        x3_idx = 5 + 8*i        # x_P_i^3
        y_idx = 6 + 8*i         # y_P_i
        y2_idx = 7 + 8*i        # y_P_i^2
        t1_idx = 8 + 8*i        # t1
        t2_idx = 9 + 8*i        # t2
        
        # Previous point indices (for i > 0)
        x_prev_idx = 3 + 8*(i-1) if i > 0 else m-2  # x_P_{i-1}
        y_prev_idx = 6 + 8*(i-1) if i > 0 else m-1  # y_P_{i-1}
        
        # 1. Constraint for k_i * (1 - k_i) = 0 (enforcing k_i is binary)
        A[constraint_idx, k_idx] = 1      # k_i
        B[constraint_idx, 0] = 1          # 1 (constant)
        B[constraint_idx, k_idx] = -1     # -k_i
        C[constraint_idx, 0] = 0          # Result should be 0
        constraint_idx += 1
        
        # 2. Constraint for x_P_i * x_P_i = x_P_i^2
        A[constraint_idx, x_idx] = 1      # x_P_i
        B[constraint_idx, x_idx] = 1      # x_P_i
        C[constraint_idx, x2_idx] = 1     # x_P_i^2
        constraint_idx += 1
        
        # 3. Constraint for x_P_i * x_P_i^2 = x_P_i^3
        A[constraint_idx, x_idx] = 1      # x_P_i
        B[constraint_idx, x2_idx] = 1     # x_P_i^2
        C[constraint_idx, x3_idx] = 1     # x_P_i^3
        constraint_idx += 1
        
        # 4. Constraint for y_P_i * y_P_i = y_P_i^2
        A[constraint_idx, y_idx] = 1      # y_P_i
        B[constraint_idx, y_idx] = 1      # y_P_i
        C[constraint_idx, y2_idx] = 1     # y_P_i^2
        constraint_idx += 1
        
        # 5. Constraint for the elliptic curve equation: y_P_i^2 = x_P_i^3 + a*x_P_i + b
        A[constraint_idx, x3_idx] = 1     # x_P_i^3
        A[constraint_idx, x_idx] = a      # a*x_P_i
        A[constraint_idx, 0] = b          # b (constant)
        B[constraint_idx, 0] = 1          # 1 (for the constant term)
        C[constraint_idx, y2_idx] = 1     # y_P_i^2
        constraint_idx += 1
        
        # 6. Constraint for computing the slope: t_1 = (δ_x*k_i + x′ - x_P_i)*(y_P_i + y_P_{i-1})
        A[constraint_idx, k_idx] = delta_x[i]  # δ_x*k_i
        A[constraint_idx, 0] = x_prime[i]      # x′
        A[constraint_idx, x_idx] = -1          # -x_P_i
        B[constraint_idx, y_idx] = 1           # y_P_i
        B[constraint_idx, y_prev_idx] = 1      # y_P_{i-1}
        C[constraint_idx, t1_idx] = 1          # t_1
        constraint_idx += 1
        
        # 7. Constraint for the other part of the slope: t_2 = (x_P_{i-1} - x_P_i)*(δ_y*k_i + y′ + y_P_i)
        A[constraint_idx, x_prev_idx] = 1      # x_P_{i-1}
        A[constraint_idx, x_idx] = -1          # -x_P_i
        B[constraint_idx, k_idx] = delta_y[i]  # δ_y*k_i
        B[constraint_idx, 0] = y_prime[i]      # y′
        B[constraint_idx, y_idx] = 1           # y_P_i
        C[constraint_idx, t2_idx] = 1          # t_2
        constraint_idx += 1
        
        # 8. Constraint for t_1 = t_2 (ensuring the slopes match)
        A[constraint_idx, t1_idx] = 1          # t_1
        B[constraint_idx, 0] = 1               # 1
        C[constraint_idx, t2_idx] = 1          # t_2
        constraint_idx += 1
    
    # Add constraint to check that the final point P_ℓ has the correct x-coordinate
    A[constraint_idx, 3 + 8*(ell-1)] = 1  # x_P_ℓ
    B[constraint_idx, 0] = 1              # 1
    C[constraint_idx, 1] = 1              # x_P (the output)
    constraint_idx += 1
    
    # Add constraint to enforce that k = sum(2^i * k_i)
    for i in range(ell):
        A[constraint_idx, 2 + 8*i] = 2^i  # 2^i * k_i
    B[constraint_idx, 0] = 1              # 1
    C[constraint_idx, 1] = 1              # k
    
    return A, B, C

# Generate a complete witness vector for the R1CS constraint system (equation 12)
def generate_witness_vector(k, H_x, G_S):
    ell = 256
    
    # Size of the witness vector
    m = 8*ell + 6
    
    # Initialize the witness vector
    z = vector(Fq, m)
    
    # First element is the constant 1
    z[0] = 1
    
    # Second element is the private key k
    z[1] = k
    
    # Third element is the x-coordinate of the final point P = k*H_x
    final_point = k * H_x
    z[2] = final_point[0]  
    
    # Define the constants c_i as per equation (9)
    c = [i+2 for i in range(ell)]
    c.append(-sum(c))  # Make the sum zero
    
    # Convert k to its binary representation (LSB first)
    k_bits = Integer(k).digits(base=2, padto=ell)
    
    # Compute all points P_i according to equation (10)
    P = []
    
    P.append(k_bits[0] * H_x + c[0] * G_S)
    
    for i in range(1, ell):

        delta_i = (2^i * k_bits[i]) * H_x + c[i] * G_S
        P.append(P[i-1] + delta_i)
    
    for i in range(ell):
        base_idx = 2 + 8*i
        
        z[base_idx] = k_bits[i]
        
        z[base_idx + 1] = P[i][0]
        
        z[base_idx + 2] = P[i][0]^2
        
        z[base_idx + 3] = P[i][0]^3
        
        z[base_idx + 4] = P[i][1]
        
        z[base_idx + 5] = P[i][1]^2
        
        # For each bit, we need to compute t_1 and t_2
        if i > 0:
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
        else:
            # For i=0, we don't have a previous point
            # Set t_1 and t_2 to 0 or some other sensible value
            z[base_idx + 6] = 0
            z[base_idx + 7] = 0
    
    # The last two elements of the witness vector are placeholders for x_P_{-1} and y_P_{-1}
    # These are used when computing t_1 and t_2 for i=0
    z[m-2] = 0  
    z[m-1] = 0  
    
    return z



    
    





  



