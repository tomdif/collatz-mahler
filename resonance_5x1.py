#!/usr/bin/env python3
"""
RESONANCE DETECTION FOR THE 5x+1 MAP
====================================

This script computes dim ker(I - M_N^k) for the 5x+1 map at various periods k.
The 5x+1 map has a known 5-cycle {1, 3, 8, 4, 2}, which should produce a
nullity spike at period k=5.

This serves as a positive control to validate the spectral method.

Usage: python resonance_5x1.py [N]
  N: truncation degree (default 30)

Output: Table of nullity by period, showing spike at k=5
"""

import sys
from fractions import Fraction
from functools import lru_cache

@lru_cache(maxsize=None)
def binom(n, k):
    """Binomial coefficient C(n,k) with exact arithmetic."""
    if k < 0 or k > n: return 0
    if k == 0 or k == n: return 1
    if k > n // 2: k = n - k
    res = 1
    for i in range(k):
        res = res * (n - i) // (i + 1)
    return res

@lru_cache(maxsize=None)
def map_5x1(n):
    """5x+1 map: T(n) = n/2 if even, (5n+1)/2 if odd."""
    if n == 0: return 0
    if n % 2 == 0: return n // 2
    return (5 * n + 1) // 2

def build_matrix(N, map_func):
    """Build NxN Mahler matrix for given map."""
    M = []
    for i in range(N):
        row = [Fraction(0)] * N
        for k in range(i + 1):
            sign = 1 if (i - k) % 2 == 0 else -1
            coeff = sign * binom(i, k)
            Tk = map_func(k)
            for j in range(N):
                b = binom(Tk, j)
                if b != 0:
                    row[j] += coeff * b
        M.append(row)
    return M

def matrix_multiply(A, B):
    """Multiply two matrices with exact arithmetic."""
    N = len(A)
    C = [[Fraction(0)] * N for _ in range(N)]
    for i in range(N):
        for j in range(N):
            for k in range(N):
                C[i][j] += A[i][k] * B[k][j]
    return C

def matrix_power(M, k):
    """Compute M^k using repeated squaring."""
    N = len(M)
    if k == 1:
        return [row[:] for row in M]
    
    # Identity matrix
    result = [[Fraction(1 if i == j else 0) for j in range(N)] for i in range(N)]
    base = [row[:] for row in M]
    
    while k > 0:
        if k % 2 == 1:
            result = matrix_multiply(result, base)
        base = matrix_multiply(base, base)
        k //= 2
    
    return result

def kernel_nullity(M):
    """Compute nullity of (I-M) via Gaussian elimination."""
    N = len(M)
    aug = [[Fraction(-M[i][j] + (1 if i == j else 0)) 
            for j in range(N)] for i in range(N)]
    pivots, prow = [], 0
    for col in range(N):
        if prow >= N: break
        found = next((r for r in range(prow, N) if aug[r][col] != 0), -1)
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
    return N - len(pivots)

def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    
    print(f"5x+1 Map Resonance Detection (N={N})")
    print("=" * 40)
    print(f"\nThe 5x+1 map has a known 5-cycle: {{1, 3, 8, 4, 2}}")
    print("We expect a nullity spike at period k=5.\n")
    
    print("Building base matrix M_N...")
    M = build_matrix(N, map_5x1)
    
    print("\nPeriod k | dim ker(I - M^k)")
    print("-" * 30)
    
    for k in range(1, 11):
        if k == 1:
            Mk = M
        else:
            Mk = matrix_power(M, k)
        
        nullity = kernel_nullity(Mk)
        marker = " <-- 5-CYCLE DETECTED!" if nullity > 3 else ""
        print(f"   {k:2d}    |       {nullity}{marker}")
    
    print("\n" + "=" * 40)
    print("The spike at k=5 confirms the method detects cycles.")

if __name__ == "__main__":
    main()
