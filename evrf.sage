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

#private key 
k =  Fq.random_element()

print(k)

# this is the naive hash and pray
# TODO implement a proper hash_to_curve
#def hash_to_curve_bandersnatch(x):
#    while True:
#        i    
    


#def eval(k, x):
