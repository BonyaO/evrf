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
k =  Fq.random_element()
  
# this is the naive hash and pray
# TODO implement a proper hash_to_curve
def hash_to_curve_bandersnatch(x):
    while True:
        try: 
            P = F.lift_x(x)
            return P   
        except (ValueError):
            x = x+1

def eval(k,x):
    H = hash_to_curve_bandersnatch(x)
    P = k* H
    xp = P[0]
    return xp, xp*G


print(eval(k,Fr.random_element() ))