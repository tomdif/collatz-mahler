# The Dipole Conjecture

**Spectral Rigidity of the Collatz Koopman Operator in the Mahler Basis**

This repository contains the code and paper for investigating the Collatz conjecture through truncated Koopman operators acting on Mahler polynomials over the 2-adic integers.

## Main Results

1. **Theorem (dim ≥ 1)**: The kernel of (I - M_N) has dimension at least 1 for all N, with the lower bound achieved by the fixed point at 0.

2. **Theorem (2/3-Rigidity)**: The alternating vector v_N = [0, 1, -1, 1, -1, ...]^T satisfies ((I - M_N) v_N)_m = 0 exactly for all m < ⌊2N/3⌋. This explains the remarkable stability of the observed kernel structure.

3. **Computational Observation**: For all N ≤ 500, dim ker(I - M_N) = 2 exactly.

4. **Dipole Conjecture**: dim ker(I - M_N) = 2 for all N.

5. **Positive Control**: The 5x+1 map (which has a known 5-cycle) shows a nullity spike from 2 to 6 at period 5, confirming the method detects cycles when they exist.

## Repository Contents

- `collatz-dipole.tex` - Full paper (LaTeX)
- `verify_dipole.py` - Main verification script for the Dipole Conjecture
- `verify_5x1_control.py` - Positive control script for 5x+1 cycle detection
- `verify_two_thirds.py` - Verification of the 2/3-Rigidity Theorem

## Quick Start

```bash
# Verify the Dipole Conjecture for N=100
python3 verify_dipole.py 100

# Run the 5x+1 positive control
python3 verify_5x1_control.py

# Verify the 2/3-Rigidity Theorem
python3 verify_two_thirds.py
```

## Requirements

- Python 3.7+
- No external dependencies (uses only `fractions`, `functools`, `hashlib` from stdlib)

## Verification Hashes (N=100)

For reproducibility, the SHA-256 hashes of the idealized kernel vectors at N=100:

```
v1 (δ₀):     33d83547e87812d859d68bc0d71edc1a970a82d714c7d170dba414ce185dd446
v2 (dipole): 95c5a2404a266b6fd867b4a054d28f7a60bc9f48785af90dda591f3137b7bff1
```

## 5x+1 Positive Control Hashes (N=30)

```
5x+1 M:    9d5cec982ae29b29d024156a437632bb
3x+1 M:    f1c9ab07df9b5ea3dfe4813aca104dd5
5x+1 M^5:  8bf59ad4c6ed21f9eb5ac6248496d922
3x+1 M^5:  853f75290b4ddb64283f0938252e2374
```

## Citation

```bibtex
@article{difiore2025dipole,
  title={The Dipole Conjecture: Spectral Rigidity of the Collatz Koopman Operator in the Mahler Basis},
  author={DiFiore, Thomas A.},
  year={2025},
  note={Preprint}
}
```

## License

MIT License - see LICENSE file.

## Author

Thomas A. DiFiore (tomdif@gmail.com)
