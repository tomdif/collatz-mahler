# Collatz-Mahler Spectral Analysis

[![Python 3.6+](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Computational verification of the **Dipole Conjecture**: the kernel of the truncated Koopman operator for the Collatz map in the Mahler basis is exactly 2-dimensional.

## The Dipole Conjecture

For the accelerated Collatz (Syracuse) map $T(n) = n/2$ if $n$ even, $(3n+1)/2$ if $n$ odd, we construct the $N \times N$ Mahler matrix $M_N$ representing the truncated Koopman operator.

**Theorem (Proved):** $\dim \ker(I - M_N) \geq 1$ for all $N$.

**Conjecture (Verified for $N \leq 500$):** $\dim \ker(I - M_N) = 2$ for all $N$.

The two kernel directions correspond to:
1. $\delta_0 = [1, 0, 0, \ldots]$ — the fixed point at 0
2. $\delta_0 - \delta_{-1} \approx [0, 1, -1, 1, -1, \ldots]$ — the fixed point at -1

## Quick Start

```bash
# Verify the Dipole Conjecture for N=100
python verify_dipole.py 100

# Positive control: detect the 5-cycle in the 5x+1 map
python resonance_5x1.py 30
```

## Files

| File | Description |
|------|-------------|
| `verify_dipole.py` | Main verification script for Collatz |
| `resonance_5x1.py` | Positive control using 5x+1 map |
| `collatz-dipole.tex` | LaTeX source of the paper |
| `collatz-dipole.pdf` | Compiled paper |

## Results

### Collatz Map (3x+1)
```
N=10:   nullity = 2 ✓
N=50:   nullity = 2 ✓
N=100:  nullity = 2 ✓
N=200:  nullity = 2 ✓
N=500:  nullity = 2 ✓
```

### 5x+1 Map (Positive Control)
The 5x+1 map has a known 5-cycle {1, 3, 8, 4, 2}. Our method detects it:

```
Period k | dim ker(I - M^k)
---------|------------------
   1     |       2
   2     |       2
   3     |       2
   4     |       2
   5     |       6  ← 5-cycle detected!
```

The nullity spike from 2 to 6 at period 5 confirms the method's sensitivity to cycles.

## Requirements

- Python 3.6+
- Standard library only (no external packages)
- Uses `fractions` module for exact rational arithmetic

## Runtime

| N | Time |
|---|------|
| 10 | ~1 second |
| 100 | ~1 minute |
| 200 | ~15 minutes |
| 500 | ~5 hours |

## Mathematical Background

The Mahler basis consists of binomial coefficient functions $\phi_n(x) = \binom{x}{n}$ which form a basis for continuous functions on $\mathbb{Z}_2$. The Koopman operator $\mathcal{K}$ acts by $(\mathcal{K}f)(x) = f(T(x))$.

The matrix entry $M[m,n]$ is computed via finite differences:
$$M[m,n] = \sum_{k=0}^{m} (-1)^{m-k} \binom{m}{k} \binom{T(k)}{n}$$

## Citation

If you use this code, please cite:

```bibtex
@article{difiore2025dipole,
  title={The Dipole Conjecture: Spectral Rigidity of the Collatz Koopman Operator in the Mahler Basis},
  author={DiFiore, Thomas A.},
  journal={Experimental Mathematics},
  year={2025},
  note={Submitted}
}
```

## License

MIT License - see [LICENSE](LICENSE) file.

## Author

Thomas A. DiFiore (tomdif@gmail.com)
