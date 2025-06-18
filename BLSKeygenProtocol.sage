# BLS12-381 key derivation from Bandersnatch using eVRF

import hashlib
import json
from typing import List, Tuple, Dict, Optional

# Load the eVRF implementations
load("evrf_common.sage")
load("basic_evrf.sage")


class BLSKeyDerivationProtocol:
    """
    Protocol for deriving multiple BLS12-381 keys from a single Bandersnatch key.
    Uses eVRF to deterministically generate BLS12-381 key pairs.
    """
    
    def __init__(self):
        """
        Initialize the protocol.
        """
        
        # BLS12-381 parameters
        self.bls_order = q  
        self.bls_g1 = G1   
        
        # Bandersnatch parameters
        self.bandersnatch_order = s  
        self.bandersnatch_g = GS     
        
    def generate_master_keypair(self) -> Tuple:
        """
        Generate a master Bandersnatch keypair.
        """
        sk, vk = keygen()
            
        return sk, vk
    
    def derive_bls_keypair(self, sk, index: int) -> Dict:
        """
        Derive a BLS12-381 keypair from master Bandersnatch key.
        """
        # Encode index as 4 bytes (big-endian)
        if index < 0 or index >= 2**32:
            raise ValueError("Index must be between 0 and 2^32-1")
            
        index_bytes = int(index).to_bytes(4, 'big')
        x = int.from_bytes(index_bytes, 'big')
        
        
        y, Y, pi = eval(sk, x)
        
       
        # Use SHA-256 to ensure uniform distribution in BLS scalar field
        y_bytes = Integer(y).to_bytes(32, 'big')
        hash_output = hashlib.sha256(y_bytes).digest()
        bls_sk = int.from_bytes(hash_output, 'big') % self.bls_order
        
        # Compute BLS public key
        bls_pk = bls_sk * self.bls_g1
        
        return {
            'bls_sk': bls_sk,
            'bls_pk': bls_pk,
            'evrf_output': y,
            'evrf_Y': Y,  
            'evrf_proof': pi,
            'index': index
        }
    
    def derive_bls_public_key(self, sk, index: int) -> Tuple:
        """
        Derive only the BLS12-381 public key (without secret key).
        """
        result = self.derive_bls_keypair(sk, index)
        return result['bls_pk'], {
            'evrf_output': result['evrf_output'],
            'evrf_Y': result['evrf_Y'],  # Include the group element Y
            'evrf_proof': result['evrf_proof'],
            'index': result['index']
        }
    
    def verify_bls_derivation(self, vk, index: int, 
                              bls_pk, evrf_output, evrf_Y, evrf_proof) -> bool:
        """
        Verify that a BLS public key was correctly derived.
        """
        # Encode index
        index_bytes = int(index).to_bytes(4, 'big')
        x = int.from_bytes(index_bytes, 'big')
        
        # Fix: Pass Y (group element) instead of evrf_output (scalar) to verify function
        if not verify(vk, x, evrf_Y, evrf_proof):
            return False
        
        # Verify that Y = evrf_output * GT2
        expected_Y = evrf_output * GT2
        if evrf_Y != expected_Y:
            return False
        
        # Recompute BLS public key from verified eVRF output
        y = evrf_output
        y_bytes = Integer(y).to_bytes(32, 'big')
        hash_output = hashlib.sha256(y_bytes).digest()
        expected_bls_sk = int.from_bytes(hash_output, 'big') % self.bls_order
        expected_bls_pk = expected_bls_sk * self.bls_g1
        
        return bls_pk == expected_bls_pk
    
    def batch_derive_bls_keys(self, sk, indices: List[int], 
                              include_private: bool = True) -> List[Dict]:
        """
        Derive multiple BLS keys in batch.
        """
        results = []
        for idx in indices:
            if include_private:
                result = self.derive_bls_keypair(sk, idx)
            else:
                bls_pk, proof = self.derive_bls_public_key(sk, idx)
                result = {
                    'bls_pk': bls_pk,
                    'index': idx,
                    **proof
                }
            results.append(result)
        return results
    
    def generate_derivation_proof(self, sk, vk, indices: List[int]) -> Dict:
        """
        Generate a proof of correct derivation for multiple BLS keys.
        """
        proof = {
            'master_vk': vk,
            'derivations': []
        }
        
        for idx in indices:
            bls_pk, derivation_data = self.derive_bls_public_key(sk, idx)
            proof['derivations'].append({
                'index': idx,
                'bls_pk': bls_pk,
                **derivation_data
            })
            
        return proof
    
    def verify_derivation_proof(self, proof: Dict) -> bool:
        """
        Verify a derivation proof for multiple BLS keys.
        """
        vk = proof['master_vk']
        
        for derivation in proof['derivations']:
            # Fix: Pass evrf_Y parameter to verify_bls_derivation
            if not self.verify_bls_derivation(
                vk,
                derivation['index'],
                derivation['bls_pk'],
                derivation['evrf_output'],
                derivation['evrf_Y'],  # Add the missing Y parameter
                derivation['evrf_proof']
            ):
                return False
                
        return True


class BLSKeyManager:
    """
    High-level manager for BLS key derivation from Bandersnatch.
    """
    
    def __init__(self):
        self.protocol = BLSKeyDerivationProtocol()
        self.master_sk = None
        self.master_vk = None
        
    def initialize_from_seed(self, seed: bytes):
        """
        Initialize master key from seed.
        """
        if len(seed) != 32:
            raise ValueError("Seed must be exactly 32 bytes")
        
        # Use seed to generate deterministic private key
        seed_int = int.from_bytes(seed, 'big')
        
        # For basic eVRF: just k
        k = seed_int % (2^(l + 1) - 1) + 1
        self.master_sk = k
        self.master_vk = k * GT1
    
    def initialize_random(self):
        """Initialize with random master key."""
        self.master_sk, self.master_vk = self.protocol.generate_master_keypair()
    
    def derive_keys(self, indices: List[int], 
                    include_private: bool = True) -> List[Dict]:
        """
        Derive BLS keys for specified indices.
        """
        if not self.master_sk:
            raise RuntimeError("Master key not initialized")
            
        return self.protocol.batch_derive_bls_keys(
            self.master_sk, indices, include_private
        )
    
    def derive_range(self, start: int, count: int, 
                     include_private: bool = True) -> List[Dict]:
        """
        Derive BLS keys for a range of indices.
        """
        indices = list(range(start, start + count))
        return self.derive_keys(indices, include_private)
    
    def export_key(self, key_data: Dict, format: str = 'decimal') -> Dict:
        """
        Export BLS key in decimal format.
        """
        export = {'index': key_data['index']}
        
        # Export public key
        pk_x = Integer(key_data['bls_pk'][0])
        export['public_key'] = pk_x
        
        # Export private key if available
        if 'bls_sk' in key_data:
            export['private_key'] = key_data['bls_sk']
                
        return export
    
    def generate_proof(self, indices: List[int]) -> Dict:
        """Generate derivation proof for specified indices."""
        if not self.master_sk:
            raise RuntimeError("Master key not initialized")
            
        return self.protocol.generate_derivation_proof(
            self.master_sk, self.master_vk, indices
        )


# Utility functions
def demonstrate_basic_usage():
    """Demonstrate basic protocol usage."""
    print("=== BLS Key Derivation Protocol Demo ===\n")
    
    # Initialize with basic eVRF
    manager = BLSKeyManager()
    
    # Generate master key from seed
    seed = b"demo-seed-for-bls-key-derivation"
    manager.initialize_from_seed(seed)
    print("✓ Master Bandersnatch key initialized\n")
    
    # print master key
    print(f"Master Bandersnatch Secret Key: {manager.master_sk}")
    print(f"Master Bandersnatch Verification Key: {manager.master_vk}\n")


    # Derive BLS keys for indices 0-1
    print("Deriving BLS keys for indices 0-1...")
    keys = manager.derive_range(start=0, count=2)
    
    for key_data in keys:
        exported = manager.export_key(key_data)
        print(f"\nIndex {exported['index']}:")
        print(f"  BLS Public Key: {str(exported['public_key'])}")
        if 'private_key' in exported:
            print(f"  BLS Private Key: {str(exported['private_key'])}")

    # Generate and verify derivation proof
    print("\n\nGenerating derivation proof for indices [0, 1]...")
    proof = manager.generate_proof([0, 1])
    
    # Verify the proof
    is_valid = manager.protocol.verify_derivation_proof(proof)
    print(f"✓ Derivation proof verified: {is_valid}")
    
    return manager


# Run demonstrations if executed directly
if __name__ == "__main__":
    print("Starting BLS Key Derivation Protocol Demo...\n")
    
    # Basic usage demo
    manager = demonstrate_basic_usage()
    
    print("\n✓ All demos completed successfully!")