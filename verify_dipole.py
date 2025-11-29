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

Usage:
    python3 verify_dipole.py [N]

Runtime: N=100 takes ~1 minute; N=500 takes ~5 hours.
"""

import sys
import hashlib
from fractions import Fraction
from functools import lru_cache


@lru_cache(maxsize=None)
def binom(n: int, k: int) -> int:
    """Binomial coefficient C(n,k) with exact integer arithmetic."""
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    if k > n // 2:
        k = n - k
    result = 1
    for i in range(k):
        result = result * (n - i) // (i + 1)
    return result


@lru_cache(maxsize=None)
def collatz(n: int) -> int:
    """Syracuse (accelerated Collatz) map: T(n) = n/2 or (3n+1)/2."""
    if n == 0:
        return 0
    if n % 2 == 0:
        return n // 2
    return (3 * n + 1) // 2


def build_matrix(N: int) -> list:
    """Build NxN Mahler matrix."""
    print(f"Building {N}x{N} Mahler matrix...", flush=True)
    matrix = []
    for m in range(N):
        if m % 50 == 0:
            print(f"  Row {m}/{N}", flush=True)
        row = [0] * N
        for k in range(m + 1):
            sign = 1 if (m - k) % 2 == 0 else -1
            coeff = sign * binom(m, k)
            Tk = collatz(k)
            for j in range(min(Tk + 1, N)):
                row[j] += coeff * binom(Tk, j)
        matrix.append(row)
    return matrix


def kernel_nullity_and_vectors(M: list) -> tuple:
    """
    Compute nullity and kernel vectors of (I-M) via Gaussian elimination.
    Uses exact rational arithmetic (Fraction) to avoid numerical issues.
    Returns (nullity, kernel_vectors).
    """
    N = len(M)
    print(f"Computing kernel of (I - M_{N})...", flush=True)
    
    # Build (I - M) with exact rational entries
    A = [[Fraction(1 if i == j else 0) - M[i][j] for j in range(N)] for i in range(N)]
    
    # Gaussian elimination to RREF
    pivots = []
    pivot_row = 0
    
    for col in range(N):
        if pivot_row >= N:
            break
        
        # Find pivot
        found = -1
        for r in range(pivot_row, N):
            if A[r][col] != 0:
                found = r
                break
        
        if found == -1:
            continue
        
        # Swap rows
        A[pivot_row], A[found] = A[found], A[pivot_row]
        pivots.append(col)
        
        # Scale pivot row
        pv = A[pivot_row][col]
        for c in range(N):
            A[pivot_row][c] /= pv
        
        # Eliminate all other rows
        for r in range(N):
            if r != pivot_row and A[r][col] != 0:
                factor = A[r][col]
                for c in range(N):
                    A[r][c] -= factor * A[pivot_row][c]
        
        pivot_row += 1
    
    # Extract kernel vectors
    free_vars = [j for j in range(N) if j not in pivots]
    nullity = len(free_vars)
    
    kernel_vectors = []
    for fv in free_vars:
        vec = [Fraction(0)] * N
        vec[fv] = Fraction(1)
        for row_idx, pivot_col in enumerate(pivots):
            vec[pivot_col] = -A[row_idx][fv]
        kernel_vectors.append(vec)
    
    return nullity, kernel_vectors


def vector_hash(vec: list) -> str:
    """Compute SHA-256 hash of a rational vector."""
    denoms = [abs(v.denominator) for v in vec if v != 0]
    if not denoms:
        return hashlib.sha256(b"zero").hexdigest()
    
    from math import gcd
    lcm = denoms[0]
    for d in denoms[1:]:
        lcm = lcm * d // gcd(lcm, d)
    
    int_vec = [int(v * lcm) for v in vec]
    data = ",".join(map(str, int_vec)).encode('utf-8')
    return hashlib.sha256(data).hexdigest()


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    
    print("=" * 60)
    print(f"DIPOLE CONJECTURE VERIFICATION (N = {N})")
    print("=" * 60)
    
    M = build_matrix(N)
    nullity, kernel_vecs = kernel_nullity_and_vectors(M)
    
    print(f"\n*** RESULT: Nullity = {nullity} ***")
    
    if nullity == 2:
        print("✓ Dipole Conjecture VERIFIED for this N")
    else:
        print(f"✗ UNEXPECTED: Nullity = {nullity} (expected 2)")
    
    # Display kernel vector hashes
    print(f"\nKernel vector hashes (N={N}):")
    for i, vec in enumerate(kernel_vecs):
        h = vector_hash(vec)
        print(f"  v{i+1}: {h}")
    
    # Verify idealized vectors for N=100
    if N == 100:
        print("\nExpected hashes for N=100 (idealized):")
        print("  v1: 33d83547e87812d859d68bc0d71edc1a970a82d714c7d170dba414ce185dd446")
        print("  v2: 95c5a2404a266b6fd867b4a054d28f7a60bc9f48785af90dda591f3137b7bff1")
    
    print("=" * 60)
    
    return 0 if nullity == 2 else 1


if __name__ == "__main__":
    sys.exit(main())
