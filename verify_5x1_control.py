#!/usr/bin/env python3
"""
POSITIVE CONTROL: 5x+1 CYCLE DETECTION
======================================

Verifies that the Mahler-Koopman spectral method correctly detects
the known 5-cycle in the 5x+1 map: {1 → 3 → 8 → 4 → 2 → 1}.

Method:
  - For fixed points: compute ker(I - M_N)
  - For p-cycles: compute ker(I - M_N^p)

Expected results:
  - 5x+1 map: nullity increases from 2 to 6 when going from M to M^5
    (4 additional vectors from the 5-cycle's complex eigenvalues)
  - 3x+1 map: nullity stays at 2 (no cycles detected)

Usage:
    python3 verify_5x1_control.py
"""

import hashlib
import sys
from fractions import Fraction
from functools import lru_cache


@lru_cache(maxsize=None)
def binom(n: int, k: int) -> int:
    """Exact binomial coefficient."""
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


def T5(n: int) -> int:
    """5x+1 Syracuse map: T(n) = n/2 if even, (5n+1)/2 if odd."""
    if n == 0:
        return 0
    return n // 2 if n % 2 == 0 else (5 * n + 1) // 2


def T3(n: int) -> int:
    """3x+1 Syracuse map (Collatz)."""
    if n == 0:
        return 0
    return n // 2 if n % 2 == 0 else (3 * n + 1) // 2


def build_matrix(N: int, T) -> list:
    """Build N×N Mahler-Koopman matrix."""
    M = [[Fraction(0)] * N for _ in range(N)]
    for m in range(N):
        for k in range(m + 1):
            sign = 1 if (m - k) % 2 == 0 else -1
            coeff = sign * binom(m, k)
            Tk = T(k)
            for j in range(min(Tk + 1, N)):
                M[m][j] += coeff * binom(Tk, j)
    return M


def mat_mult(A: list, B: list) -> list:
    """Matrix multiplication with exact arithmetic."""
    N = len(A)
    C = [[Fraction(0)] * N for _ in range(N)]
    for i in range(N):
        for j in range(N):
            for k in range(N):
                C[i][j] += A[i][k] * B[k][j]
    return C


def mat_power(M: list, p: int) -> list:
    """Compute M^p via binary exponentiation."""
    N = len(M)
    if p == 1:
        return [row[:] for row in M]
    result = [[Fraction(1 if i == j else 0) for j in range(N)] for i in range(N)]
    base = [row[:] for row in M]
    while p > 0:
        if p % 2 == 1:
            result = mat_mult(result, base)
        base = mat_mult(base, base)
        p //= 2
    return result


def kernel_nullity(M: list, power: int = 1) -> int:
    """Compute nullity of (I - M^power)."""
    N = len(M)
    Mp = mat_power(M, power) if power > 1 else M
    
    # Build I - M^p
    A = [[Fraction(1 if i == j else 0) - Mp[i][j] for j in range(N)] for i in range(N)]
    
    # Gaussian elimination
    pivot_row = 0
    for col in range(N):
        found = -1
        for r in range(pivot_row, N):
            if A[r][col] != 0:
                found = r
                break
        if found == -1:
            continue
        A[pivot_row], A[found] = A[found], A[pivot_row]
        pv = A[pivot_row][col]
        for c in range(N):
            A[pivot_row][c] /= pv
        for r in range(N):
            if r != pivot_row and A[r][col] != 0:
                factor = A[r][col]
                for c in range(N):
                    A[r][c] -= factor * A[pivot_row][c]
        pivot_row += 1
    
    return N - pivot_row


def vector_hash(vec: list) -> str:
    """SHA-256 hash of rational vector."""
    denoms = [abs(v.denominator) for v in vec if v != 0]
    if not denoms:
        return hashlib.sha256(b"zero").hexdigest()
    from math import gcd
    lcm = denoms[0]
    for d in denoms[1:]:
        lcm = lcm * d // gcd(lcm, d)
    data = ",".join(str(int(v * lcm)) for v in vec).encode()
    return hashlib.sha256(data).hexdigest()


def matrix_hash(M: list) -> str:
    """Hash of matrix entries."""
    flat = [M[i][j] for i in range(len(M)) for j in range(len(M))]
    return vector_hash(flat)[:32]


def main():
    N = 30  # Matrix power computation is O(N^3 log p), keep moderate
    
    print("=" * 70)
    print("POSITIVE CONTROL: 5x+1 CYCLE DETECTION")
    print("=" * 70)
    
    # Show the 5-cycle
    print("\n5x+1 map has known 5-cycle: 1 → 3 → 8 → 4 → 2 → 1")
    
    # Build matrices
    print(f"\nBuilding matrices (N = {N})...")
    M5 = build_matrix(N, T5)
    M3 = build_matrix(N, T3)
    
    # Compute nullities
    print("Computing kernel dimensions...")
    
    results = {}
    for name, M in [("5x+1", M5), ("3x+1", M3)]:
        results[name] = {}
        for p in [1, 5]:
            print(f"  {name} ker(I - M^{p})...", end=" ", flush=True)
            null = kernel_nullity(M, p)
            results[name][p] = null
            print(f"nullity = {null}")
    
    # Summary
    print("\n" + "=" * 70)
    print("RESULTS")
    print("=" * 70)
    print(f"\n{'Map':<8} | {'ker(I-M)':<12} | {'ker(I-M^5)':<12} | {'Δ':<6}")
    print("-" * 45)
    for name in ["5x+1", "3x+1"]:
        n1, n5 = results[name][1], results[name][5]
        delta = n5 - n1
        print(f"{name:<8} | {n1:<12} | {n5:<12} | {delta:<6}")
    
    print("\n" + "=" * 70)
    print("INTERPRETATION")
    print("=" * 70)
    
    delta_5x1 = results["5x+1"][5] - results["5x+1"][1]
    delta_3x1 = results["3x+1"][5] - results["3x+1"][1]
    
    print(f"""
5x+1: Nullity increases by {delta_5x1} when computing M^5
      → The 5-cycle {{1,2,3,4,8}} creates 4 additional eigenvectors
      → These correspond to 5th roots of unity (e^{{2πik/5}})

3x+1: Nullity stays constant (Δ = {delta_3x1})
      → No 5-cycles detected
      → Consistent with Collatz conjecture
""")
    
    # Verification hashes
    print("=" * 70)
    print("SHA-256 VERIFICATION HASHES")
    print("=" * 70)
    
    print(f"\nMatrix hashes (N={N}):")
    print(f"  5x+1 M:    {matrix_hash(M5)}")
    print(f"  3x+1 M:    {matrix_hash(M3)}")
    
    M5_5 = mat_power(M5, 5)
    M3_5 = mat_power(M3, 5)
    print(f"  5x+1 M^5:  {matrix_hash(M5_5)}")
    print(f"  3x+1 M^5:  {matrix_hash(M3_5)}")
    
    print("\n" + "=" * 70)
    success = delta_5x1 > 0 and delta_3x1 == 0
    if success:
        print("✓ POSITIVE CONTROL VALIDATED")
        print("  The method correctly detects cycles when they exist.")
    else:
        print("✗ UNEXPECTED RESULT")
    print("=" * 70)
    
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
