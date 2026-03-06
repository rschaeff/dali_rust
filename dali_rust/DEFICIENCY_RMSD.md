# Deficiency Report: RMSD Always Returns 9.9

**Reported**: 2026-03-06
**Resolved**: 2026-03-06
**Severity**: Moderate — affects result interpretation, not search accuracy
**Component**: `dali-core/src/dp/mod.rs`, line 166

## Summary

The RMSD field in all DALI search results is hardcoded to `9.9` and never actually computed. This was discovered during batch structural searches of 38 archaeal dark matter proteins against the PDB chain library (27,289 targets). Every hit, regardless of Z-score (ranging from 2.3 to 26.6), reported RMSD = 9.9.

## Evidence

All 32 hits across 38 completed search jobs report identical RMSD:

```
T1_C0141  4av3A  Z=26.6  RMSD=9.9  nblock=6   (strong hit, should have low RMSD)
T1_C0007  8qadA  Z=11.0  RMSD=9.9  nblock=8
T1_C0136  2eo0A  Z=7.5   RMSD=9.9  nblock=10
T1_C0005  2zbkD  Z=3.7   RMSD=9.9  nblock=5   (weak hit, RMSD should be higher)
```

A Z-score of 26.6 with 6 aligned blocks should produce an RMSD in the 1-3 Å range, not 9.9.

## Root Cause

In `dali-core/src/dp/mod.rs`, the `DccpEntry` is constructed with a placeholder:

```rust
// dp/mod.rs, line 160-168
if zmax >= zcut {
    results.push(DccpEntry {
        cd1: cd1.clone(),
        cd2: cd2.clone(),
        score: x1,
        zscore: zmax,
        rmsd: 9.9,          // <-- hardcoded placeholder, never computed
        blocks: blocks.clone(),
    });
}
```

## DaliLite v5 Reference

The original Fortran source (`comparemodules.f`, lines 1138-1139) has the same gap — `getrmsd()` is commented out at the DCCP stage:

```fortran
        ! rmsd is filled by fssp program
!       call getrmsd(rmsd,nblock,l1,r1,l2,r2)
```

However, the Fortran code does contain a working `getrmsd` subroutine (lines 1160-1199) that:
1. Extracts aligned Cα coordinates from both proteins using the block ranges
2. Calls `u3b()` (Kabsch/SVD superposition) with uniform weights
3. Computes `rmsd = sqrt(ssq / n)` where `ssq` is the sum of squared distances after optimal superposition
4. Returns `-9.9` on error (mismatched alignment lengths)

## Proposed Fix

Implement RMSD computation in the Rust `dowork_dp` function (or a new helper) following the Fortran `getrmsd` logic:

1. For each hit that passes the Z-score cutoff, collect Cα coordinates from both proteins using the alignment block ranges (`l1[i]..r1[i]` → protein 1, `l2[i]..r2[i]` → protein 2)
2. Verify equal number of aligned residues from both sides
3. Compute optimal superposition using SVD (Kabsch algorithm) — the rotation/translation matrices are already being computed elsewhere for the `rotation` and `translation` fields
4. Return `sqrt(sum_squared_distances / n_aligned)`

Since the `DccpEntry` already carries the alignment blocks and both proteins' coordinates are available in the store, this should be straightforward. The existing `rotation`/`translation` computation could potentially be reused or extended to also yield the RMSD.

## Resolution

Added `calc_rmsd()` to `dali-core/src/numerics/kabsch.rs` — extracts aligned CA
coordinates from alignment blocks, runs u3b (Kabsch SVD) with uniform weights,
returns `sqrt(ssq / n)`. The existing `extract_aligned_coords` helper (refactored
from `compute_transform`) avoids code duplication.

`run_dp()` now calls `calc_rmsd()` for each hit that passes the Z-score cutoff.

Verified on reference pairs:
```
101mA vs 1a00A: z=20.6, rmsd=1.52 Å (strong globin-globin)
101mA vs 1binA: z=14.5, rmsd=2.57 Å (moderate cross-fold)
1a87A vs 1allA: z= 7.3, rmsd=2.88 Å (weak match)
```

## Impact (prior to fix)

- **Search results**: Z-scores and hit detection were unaffected (RMSD is not used in scoring)
- **Result display**: RMSD was shown to users in the rustdali_server web UI and API responses
- **Downstream analysis**: Any analysis relying on RMSD values from dali_rust results got meaningless data
- **Comparison with DaliLite**: Results could not be directly compared on the RMSD dimension
