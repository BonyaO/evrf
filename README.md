# evrf
## Overview
This project implements a DDH-based Exponent Verifiable Random Function (EVRF) using Bulletproofs technology. The implementation is done in SageMath and utilizes both BLS12-381 and Bandersnatch curves.

## Components
- **Basic EVRF**: Core implementation of the EVRF functionality
- **Full EVRF**: Extended implementation with additional features
- **BLSKeygenProtocol**: Module for generating multiple public keys from a single Bandersnatch key
- **Common**: Shared utilities and Bulletproof implementation

## Usage
1. Ensure SageMath is installed on your system
2. Clone the repository
3. Run individual modules using:
```bash
sage <module_name>.sage
```

## References
The implementation is based on the following research papers:

1. ["Bulletproofs for R1CS: Bridging the Completeness-Soundness Gap and a ZK Extension"](https://eprint.iacr.org/2025/327.pdf)
    - Author: Gil Segev

2. ["Exponent-VRFs and Their Applications"](https://link.springer.com/chapter/10.1007/978-3-031-91098-2_8)
    - Authors: Dan Boneh, Iftach Haitner, Yehuda Lindell, Gil Segev

