#!/usr/bin/env python3
"""
COLLATZ MAHLER MATRIX VERIFICATION
==================================

Verifies the Dipole Conjecture: dim ker(I - M_N) = 2 for all N.

The matrix M_N represents the truncated Koopman operator for the 
accelerated Collatz map in the Mahler basis. Entry M[m,n] is computed
via the finite difference formula:
    M[m,n] = sum_{k=0}^{m} (-1)^{m-k} C(m,k) C(T(k), n)
where T is the Syracuse map and C(a,b) is the binomial coefficient.

Usage: python verify_dipole.py [N]
  N: truncation degree (default 100)

Runtime estimates:
  N=100:  ~1 minute
  N=200:  ~15 minutes  
  N=500:  ~5 hours

Output: Nullity of (I - M_N), which should equal 2.
"""

import sys
from fractions import Fraction
from functools import lru_cache
import multiprocessing
import time

@lru_cache(maxsize=None)
def binom(n, k):
    """Binomial coefficient C(n,k) with exact integer arithmetic."""
    if k < 0 or k > n: return 0
    if k == 0 or k == n: return 1
    if k > n // 2: k = n - k  # Use symmetry for efficiency
    res = 1
    for i in range(k):
        res = res * (n - i) // (i + 1)
    return res

@lru_cache(maxsize=None)
def collatz(n):
    """Syracuse (accelerated Collatz) map: T(n) = n/2 or (3n+1)/2."""
    if n == 0: return 0
    if n % 2 == 0: return n // 2
    return (3 * n + 1) // 2

def compute_row(args):
    """Compute row m of Mahler matrix M_N using finite differences."""
    i, N = args
    row = [0] * N
    # Finite difference coefficients: (-1)^{m-k} * C(m,k)
    coeffs = [(1 if (i-k) % 2 == 0 else -1) * binom(i, k) 
              for k in range(i + 1)]
    Tk = [collatz(k) for k in range(i + 1)]
    for j in range(N):
        row[j] = sum(c * binom(Tk[k], j) 
                     for k, c in enumerate(coeffs) 
                     if binom(Tk[k], j) != 0)
    return i, row

def build_matrix(N):
    """Build NxN Mahler matrix using parallel computation."""
    matrix = [None] * N
    with multiprocessing.Pool() as pool:
        for i, row in pool.imap_unordered(
            compute_row, [(i, N) for i in range(N)]):
            matrix[i] = row
    return matrix

def kernel_nullity(M):
    """Compute nullity of (I-M) via Gaussian elimination.
    
    Uses exact rational arithmetic (Fraction) to avoid 
    numerical precision issues. Returns N - rank(I-M).
    """
    N = len(M)
    # Build (I - M) with exact rational entries
    aug = [[Fraction(-M[i][j] + (1 if i == j else 0)) 
            for j in range(N)] for i in range(N)]
    pivots, prow = [], 0
    # Standard row reduction
    for col in range(N):
        if prow >= N: break
        found = next((r for r in range(prow, N) 
                     if aug[r][col] != 0), -1)
        if found == -1: continue
        aug[prow], aug[found] = aug[found], aug[prow]
        pivots.append(col)
        pv = aug[prow][col]
        for c in range(col, N): aug[prow][c] /= pv
        for r in range(prow + 1, N):
            if aug[r][col] != 0:
                f = aug[r][col]
                for c in range(col, N): aug[r][c] -= f * aug[prow][c]
        prow += 1
    return N - len(pivots)  # nullity = N - rank

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    
    print(f"Collatz Dipole Conjecture Verification")
    print(f"=" * 40)
    print(f"Truncation degree N = {N}")
    print()
    
    print(f"Building Mahler matrix M_{N}...")
    start = time.time()
    M = build_matrix(N)
    build_time = time.time() - start
    print(f"  Matrix built in {build_time:.1f} seconds")
    
    print(f"Computing nullity of (I - M_{N})...")
    start = time.time()
    nullity = kernel_nullity(M)
    kernel_time = time.time() - start
    print(f"  Kernel computed in {kernel_time:.1f} seconds")
    
    print()
    print(f"RESULT: N={N}, nullity = {nullity}")
    
    if nullity == 2:
        print("✓ Dipole Conjecture VERIFIED for this N")
    else:
        print("✗ Dipole Conjecture VIOLATED!")
    
    assert nullity == 2, f"Dipole Conjecture violated at N={N}!"
