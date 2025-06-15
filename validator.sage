load("evrf_common.sage")
load("basic_evrf.sage")

class BLSKeyGeneratorFromBandersnatch:
    def __init__(self, sk_bandersnatch):
        """
        Takes a Bandersnatch private key (scalar in Z_bandersnatch)
        """
        self.sk = sk_bandersnatch
        self.pk = self.sk * G  # G is the Bandersnatch group generator

    def derive_validator_output(self, index):
        """
        Use the DDH-based EVRF to derive y and π
        msg = index encoded as 4-byte big-endian
        """
        msg = int(index).to_bytes(4, 'big')
        y, pi = eval(self.sk, msg)
        return msg, y, pi

    def verify_vrf_output(self, index, y, pi):
        msg = int(index).to_bytes(4, 'big')
        return verify(self.pk, msg, y, pi)

    def derive_bls_sk(self, index):
        """
        Maps the EVRF output y to a scalar in Z_q (BLS12-381 secret key space)
        """
        _, y, _ = self.derive_validator_output(index)
        y_bytes = Integer(y).to_bytes(32, 'big')
        bls_sk = Integer.from_bytes(hashlib.sha256(y_bytes).digest()) % q  # q is the BLS scalar field order
        return bls_sk

    def derive_bls_pk(self, index):
        """
        Derives the BLS12-381 public key from the Bandersnatch private key
        index: integer ID for validator
        Output: G1 element = bls_sk * G1 (BLS generator)
        """
        bls_sk = self.derive_bls_sk(index)
        bls_pk = bls_sk * G1  # G1 is the BLS12-381 G1 generator from evrf_common.sage
        return bls_pk


def test_bls_keygen_from_bandersnatch():
    """
    Test the BLS key generation from Bandersnatch private key
    """
    sk_bandersnatch = Integer(123456789)  # Example Bandersnatch private key
    index = 42  # Example validator index

    keygen = BLSKeyGeneratorFromBandersnatch(sk_bandersnatch)
    
    # Derive BLS secret key
    bls_sk = keygen.derive_bls_sk(index)
    print(f"BLS Secret Key: {bls_sk}")

    # Derive BLS public key
    bls_pk = keygen.derive_bls_pk(index)
    print(f"BLS Public Key: {bls_pk}")