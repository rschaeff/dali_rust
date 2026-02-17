# Benchmark Results: Rust dali-core vs Fortran DaliLite.v5

**Date:** 2026-02-17
**Corpus:** 5-structure globin set (101mA, 1a00A, 1a87A, 1allA, 1binA)
**Hardware:** leda (shared server, Linux 5.15.0-88-generic, x86-64)
**Rust:** opt-level=3, lto="fat", codegen-units=1
**Fortran:** gfortran -O3 (serialcompare binary, chained via Perl as dali.pl does)

## Head-to-Head: Single-Pair Timing

Rust process vs Fortran process, best of 3 runs, same .dat files.

### WOLF (SSE hashing) — kernel only

| Pair | Nres | Rust (ms) | Fortran (ms) | Ratio | Winner |
|------|------|-----------|--------------|-------|--------|
| 101mA/1a00A | 295 | 1.2 | 2.0 | 0.59x | Rust |
| 101mA/1binA | 297 | 3.3 | 3.0 | 1.07x | Fortran |
| 1a87A/1allA | 457 | 0.5 | 1.8 | 0.28x | Rust |
| 1a87A/101mA | 451 | 0.4 | 1.7 | 0.21x | Rust |
| 1a87A/1binA | 440 | 0.4 | 2.0 | 0.20x | Rust |
| 101mA/1allA | 314 | 1.7 | 2.5 | 0.69x | Rust |
| **TOTAL** | | **7.5** | **13.1** | **0.57x** | **Rust** |

Rust WOLF is 1.8x faster overall. Fortran times cluster at ~2ms, suggesting
process startup overhead dominates. Rust times vary 0.4-3.3ms depending on
actual computational work (spatial hash collisions).

### PARSI (branch-and-bound search) — kernel only

| Pair | Nres | Rust (ms) | Fortran (ms) | Ratio | Winner |
|------|------|-----------|--------------|-------|--------|
| 101mA/1a00A | 295 | 51.2 | 20.7 | 2.47x | Fortran |
| 101mA/1binA | 297 | 47.8 | 18.7 | 2.55x | Fortran |
| 1a87A/1allA | 457 | 125.2 | 57.8 | 2.17x | Fortran |
| 1a87A/101mA | 451 | 121.5 | 51.6 | 2.36x | Fortran |
| 1a87A/1binA | 440 | 93.2 | 42.7 | 2.18x | Fortran |
| 101mA/1allA | 314 | 44.3 | 18.4 | 2.41x | Fortran |
| **TOTAL** | | **483.4** | **210.0** | **2.30x** | **Fortran** |

Fortran PARSI is 2.3x faster. PARSI is compute-bound (integer scoring in
tight loops), where Fortran's advantages dominate: polynomial eigenvalue SVD,
integer*2 distances, no bounds checking, COMMON blocks (zero allocation).

### Full Pipeline (both directions, all stages)

Fortran pipeline chains stages exactly as `dali.pl` does via `mpidali.pm`:
- Wolf path:  `serialcompare WOLF` → `serialcompare DP` → `dccp2dalicon.pl` → `serialcompare DALICON T` → `serialcompare DP`
- Parsi path: `serialcompare PARSI` → `serialcompare FILTER95` → `sort|uniq|pipe96-free.pl` → `serialcompare DALICON T` → `serialcompare DP`

Both paths run in both directions (forward + reverse), matching Rust `compare_pair`.

| Pair | Nres | Rust (ms) | Fortran (ms) | Ratio | Winner |
|------|------|-----------|--------------|-------|--------|
| 101mA/1a00A | 295 | 657 | 286 | 2.30x | Fortran |
| 101mA/1binA | 297 | 507 | 148 | 3.44x | Fortran |
| 1a87A/1allA | 457 | 1676 | 289 | 5.80x | Fortran |
| 1a87A/101mA | 451 | 1527 | 230 | 6.63x | Fortran |
| 1a87A/1binA | 440 | 1468 | 176 | 8.35x | Fortran |
| 101mA/1allA | 314 | 605 | 194 | 3.11x | Fortran |
| **TOTAL** | | **6440** | **1323** | **4.87x** | **Fortran** |

Fortran full pipeline is **4.9x faster** despite the overhead of 8-10
subprocess spawns, file I/O, and Perl format conversion per pair. Note that the
Fortran pipeline includes these subprocess/Perl costs — they are not subtracted.

The ratio ranges from 2.3x (small same-fold pair) to 8.4x (medium cross-fold
pair), correlating with how many DALICON iterations run.

## Rust Kernel Breakdown (Criterion)

Microbenchmarks on individual kernels, tested on 3 protein pairs spanning
small (141+154 res), medium (297+160 res), and large (456+378 res).

### Numerics

| Kernel | Small | Medium | Large | Scaling |
|--------|-------|--------|-------|---------|
| u3b (Kabsch SVD) | 3.3 us | 3.8 us | 7.7 us | ~O(n) |
| dpgetdist (dist matrix) | 72 us (154r) | 267 us (297r) | 632 us (456r) | O(n^2) |
| fitz (iterative superpos) | 185 us | 405 us | 5.2 ms | O(n * iters) |

### Module-Level

| Module | Small | Medium | Large | Notes |
|--------|-------|--------|-------|-------|
| wolf_compare | 884 us | 105 us | 14.4 ms | Varies by hash collisions |
| run_dp | 162 us | - | 388 us | Fast (post-WOLF) |
| run_dalicon | 65.9 ms | - | - | GA refinement, many fitz calls |
| dowork_parsi | 51.6 ms | 127 ms | 288 ms | Dominant cost |
| run_filter95 | 963 us | 12.2 ms | 100 ms | Fitz-heavy for large inputs |

### Pipeline-Level

| Path | Small | Medium | Large |
|------|-------|--------|-------|
| Wolf path (WOLF→DP→DALICON→DP) | 67.6 ms | 432 us | 15.5 ms |
| Parsi path (PARSI→FILTER95→DALICON→DP) | 118 ms | 248 ms | 424 ms |
| **Full compare_pair** | **665 ms** | **1.70 s** | **7.57 s** |

## Analysis

### Cost Attribution

For a typical small-globin pair (101mA/1a00A), Rust takes ~660ms. Where does it go?

| Stage | Time | Fraction |
|-------|------|----------|
| WOLF (2 directions) | ~2 ms | <1% |
| PARSI (2 directions) | ~100 ms | 15% |
| FILTER95 (2 directions) | ~2 ms | <1% |
| DP (8 calls total) | ~1 ms | <1% |
| **DALICON (4 calls total)** | **~550 ms** | **83%** |

DALICON dominates. Each call runs the genetic algorithm with iterative fitz
superposition. At ~66ms per call for small proteins (more for larger ones),
4 calls (forward WOLF, forward PARSI, reverse WOLF, reverse PARSI) account
for the vast majority of runtime.

### Where Rust Loses Time

1. **DALICON is 83% of pipeline cost.** The genetic algorithm runs many
   fitz+u3b iterations. Each fitz call is 185us-5.2ms depending on protein
   size, and DALICON calls fitz dozens of times per invocation.

2. **PARSI is 2.3x slower than Fortran kernel-to-kernel.** The inner scoring
   loop uses bounds-checked array access and i32 arithmetic where Fortran uses
   integer*2. Fortran's integer scoring table is a single COMMON block array;
   Rust allocates Vec per function call.

3. **Pipeline gap widens with protein size.** 2.3x for 295-res pairs but
   8.4x for 440-res pairs, because DALICON cost scales super-linearly with
   protein size (more GA iterations, larger fitz problems).

### Where Rust is Competitive

1. **WOLF is 1.8x faster than Fortran** at the kernel level. The spatial hash
   grid is a natural fit for Rust's data structures.

2. **u3b (Kabsch):** 3.3 us for 150 residues is fast. Not the bottleneck.

3. **dpgetdist:** Distance matrix computation at 72-632 us is reasonable for
   O(n^2) work.

### Optimization Opportunities (ordered by expected impact)

1. **DALICON profiling** — Since DALICON is 83% of cost, even a 2x speedup
   there would cut total pipeline time nearly in half. Profile to find whether
   the bottleneck is fitz iterations, u3b calls, or the GA loop itself.

2. **Rayon parallelism** — Forward and reverse paths can run in parallel for
   a potential 2x improvement. Fortran's serial architecture (COMMON blocks,
   single process) cannot do this.

3. **PARSI inner loop** — `unsafe { get_unchecked() }` on the hot paths in
   `segsegscore` (5 nested loops with bounds checking) could close the 2.3x
   kernel gap.

4. **Allocation elimination** — `totscore()` allocates a Vec per call.
   `transrotate()` clones columns via `to_owned()`. Pre-allocating scratch
   buffers would help.

## Summary

| Level | Ratio | Notes |
|-------|-------|-------|
| WOLF kernel | **Rust 1.8x faster** | Spatial hashing favors Rust |
| PARSI kernel | **Fortran 2.3x faster** | Integer scoring, no bounds checks |
| Full pipeline | **Fortran 4.9x faster** | DALICON (GA+fitz) dominates |

The Fortran pipeline's 4.9x advantage is real but includes 8-10 subprocess
spawns and Perl text processing per pair. In a production setting where the
Fortran overhead is amortized across batch comparisons, the gap may be larger.
Conversely, Rayon parallelism (forward+reverse in parallel) could halve the
Rust time, closing to ~2.5x.

## Fortran .dat Compatibility Note

The 18-structure ECOD corpus .dat files cause Fortran `serialcompare` to segfault
in `setupprotein_`. All benchmarks use the original 5-structure corpus whose .dat
files were generated by Fortran's own pipeline.
