#!/usr/bin/env python3
"""
VERIFICATION OF THE 2/3-RIGIDITY THEOREM
========================================

This script verifies the 2/3-Rigidity Theorem:

    THEOREM: Let v_N = [0, 1, -1, 1, -1, ...]^T. Then
             ((I - M_N) v_N)_m = 0 for all m < floor(2N/3).

The theorem is PROVED in the paper using the binomial theorem.
This script provides computational verification.

The threshold 2N/3 arises because:
  - For odd k, T(k) = (3k+1)/2 ≈ (3/2)k
  - The inner sum equals 1 when T(k) < N
  - This fails when k ≈ 2N/3, so T(k) ≈ N

Usage:
    python3 verify_two_thirds.py
"""

import sys
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


def T(n: int) -> int:
    """Syracuse map."""
    if n == 0:
        return 0
    if n % 2 == 0:
        return n // 2
    return (3 * n + 1) // 2


def compute_error_vector(N: int) -> list:
    """
    Compute (I - M_N) @ v_alt where v_alt = [0, 1, -1, 1, -1, ...].
    Returns the error vector.
    """
    v_alt = [0] + [(-1) ** (j + 1) for j in range(1, N)]
    error = []
    
    for m in range(N):
        # Compute [M_N @ v_alt]_m
        Mv_m = 0
        for k in range(m + 1):
            sign = 1 if (m - k) % 2 == 0 else -1
            coeff = sign * binom(m, k)
            Tk = T(k)
            # Inner product with v_alt
            inner = 0
            for j in range(1, min(Tk + 1, N)):
                inner += binom(Tk, j) * v_alt[j]
            Mv_m += coeff * inner
        
        # Error = v_alt[m] - [M_N @ v_alt]_m
        error.append(v_alt[m] - Mv_m)
    
    return error


def first_nonzero_index(error: list) -> int:
    """Find the first non-zero entry in the error vector."""
    for i, val in enumerate(error):
        if val != 0:
            return i
    return len(error)


def main():
    print("=" * 70)
    print("VERIFICATION OF THE 2/3-RIGIDITY THEOREM")
    print("=" * 70)
    print()
    print("THEOREM: ((I - M_N) v_alt)_m = 0 for all m < floor(2N/3)")
    print()
    print(f"{'N':>5} | {'floor(2N/3)':>12} | {'First ≠ 0':>12} | {'Verified':>10}")
    print("-" * 50)
    
    all_verified = True
    
    for N in [10, 20, 30, 50, 75, 100, 150, 200, 300]:
        error = compute_error_vector(N)
        first_nz = first_nonzero_index(error)
        threshold = (2 * N) // 3
        
        # The theorem says error[m] = 0 for m < threshold
        # So first_nz should be >= threshold
        verified = first_nz >= threshold
        if not verified:
            all_verified = False
        
        status = "✓" if verified else "✗"
        print(f"{N:>5} | {threshold:>12} | {first_nz:>12} | {status:>10}")
    
    print("-" * 50)
    print()
    
    if all_verified:
        print("✓ 2/3-RIGIDITY THEOREM VERIFIED FOR ALL TESTED N")
    else:
        print("✗ VERIFICATION FAILED")
    
    # Show detailed analysis for N=100
    print()
    print("=" * 70)
    print("DETAILED ANALYSIS FOR N = 100")
    print("=" * 70)
    
    N = 100
    error = compute_error_vector(N)
    threshold = (2 * N) // 3
    first_nz = first_nonzero_index(error)
    
    print(f"\nTheoretical threshold: floor(2*100/3) = {threshold}")
    print(f"First non-zero index: {first_nz}")
    print(f"First non-zero value: {error[first_nz]}")
    print()
    
    print("Error vector around the threshold:")
    for m in range(max(0, first_nz - 3), min(N, first_nz + 5)):
        marker = " <-- first non-zero" if m == first_nz else ""
        print(f"  error[{m}] = {error[m]}{marker}")
    
    print()
    print("=" * 70)
    print("EXPLANATION")
    print("=" * 70)
    print("""
The 2/3 threshold arises from the expansion factor of the Collatz map:

  - For odd k: T(k) = (3k+1)/2 ≈ (3/2)k
  - The inner sum S_k = Σ_j C(T(k),j) * (-1)^{j+1} equals 1 when T(k) < N
  - When T(k) ≥ N, the sum is truncated and S_k ≠ 1
  - This truncation first occurs when k ≈ 2N/3

Therefore: ((I - M_N) v_alt)_m = 0 exactly for m < 2N/3.

This is proved rigorously in the paper using the binomial theorem.
""")
    
    return 0 if all_verified else 1


if __name__ == "__main__":
    sys.exit(main())
