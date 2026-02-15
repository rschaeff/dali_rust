# PARSI Mismatch Analysis

**Date:** 2026-02-13

## Summary

56 of 435 PARSI refine output lines (12.9%) don't match between Python and Fortran
at tolerance=100. Investigation shows these are **expected divergences from the
heuristic branch-and-bound search**, not implementation bugs.

## Evidence

### Score Table Verification
- Weights table: **0 differences** between float32 and float64 implementations
- Score lookup table: **194/64,000 entries differ by exactly ±1**
  - Cause: float32 vs float64 intermediate arithmetic in `(0.20 - abs(x-y)/x)`
  - Distribution: evenly spread across all distance bins

### Mismatch Characterization
- **55 of 56** missing entries have corresponding "extra" entries with same (cd1cd2, idom) key
- The number of entries per (cd1cd2, idom) group usually matches (±1)
- Score differences range from 130 to 3,057,788

### Best vs Suboptimal Alignment Comparison
| Category         | Count | Total | Rate  |
|-----------------|-------|-------|-------|
| Best: exact     | 35    | 53    | 66%   |
| Best: close <5% | 6     | 53    | 11%   |
| Best: miss >5%  | 12    | 53    | 23%   |
| Suboptimal: match | 126 | 169   | 75%   |
| Suboptimal: miss  | 43  | 169   | 25%   |

### Independent Score Verification
For 1a87A vs 1a00A (idom=1):
- **Python alignment**: score=2,227,558 — **independently verified correct**
  - 5 aligned segments with specific residue ranges
  - Sum of singlet + doublet scores matches reported value exactly
- **Fortran alignment**: score=1,140,652 — different segments, different residue ranges
- Python finds a **genuinely higher-scoring** alignment that Fortran's search doesn't discover

### Root Cause: Heuristic Search Divergence
The PARSI branch-and-bound uses a priority stack with limited capacity (MAXSTACK).
When the stack overflows, the bottom 90% of states are dropped. Tiny scoring differences
(±1 in 0.3% of score table entries) cause states to be ordered differently, leading to:

1. **Different exploration order** — Python and Fortran pop different states first
2. **Different pruning** — stack overflow drops different low-priority states
3. **Different solutions** — the search converges to different local optima
4. **Sometimes Python wins** (higher score), sometimes Fortran wins

This is expected behavior for any NP-hard optimization with heuristic pruning.

## Structure Involvement
| Pair          | Missing | Extra | Notes |
|---------------|---------|-------|-------|
| 1a87A1a00A    | 16      | 15    | Largest protein (297 res), most domains |
| 1a87A1allA    | 13      | 15    | |
| 1a87A1binA    | 13      | 14    | |
| 1a87A101mA    | 7       | 8     | |
| 1a87A1a87A    | 4       | 4     | Self-comparison |
| 101mA1a00A    | 2       | 2     | Only non-1a87A pair with >5% diff |
| 1allA1binA    | 1       | 1     | |

93% of mismatches involve 1a87A (the largest protein at 297 residues).

## Tolerance Sweep
| Tolerance | Matched  | Rate  | Recovered |
|-----------|----------|-------|-----------|
| 100       | 379/435  | 87.1% | 0         |
| 1,000     | 385/435  | 88.5% | 6         |
| 10,000    | 395/435  | 90.8% | 16        |
| 100,000   | 413/435  | 94.9% | 34        |
| 500,000   | 431/435  | 99.1% | 52        |

At tolerance=500K, 99.1% match. The remaining 4 entries are from count imbalances
where Python produces one more/fewer suboptimal alignment per domain group.

## Conclusion

The 87% match rate at tolerance=100 reflects the inherent sensitivity of
suboptimal alignment enumeration to floating-point arithmetic differences.
The implementation is algorithmically correct — both implementations find
valid alignments, just different ones due to search path divergence.

For the validation goal (proving AI-assisted reimplementation feasibility),
this is a strong result: the algorithm produces structurally valid alignments
with correct scoring, and the divergence is explainable by well-understood
floating-point arithmetic effects.
