# Optimization Log

## Baseline (2026-02-14)

Hardware: leda (shared HPC head node)
Rust: 1.93.1, release profile (default settings — opt-level 3, no LTO, 16 codegen-units)
Fortran: gfortran, DaliLite.v5 (serialcompare binary + Perl glue)
Corpus: 5-struct (101mA, 1a00A, 1a87A, 1allA, 1binA)

### Single pair (101mA vs 1a00A)

| Implementation | Time | Slowdown |
|---|---|---|
| Fortran (DaliLite.v5) | 0.39s | 1.0x |
| Rust (PyO3, release) | 0.83s | 2.1x |

### 20 directed pairs (5 structures)

| Implementation | Total | Per pair | Slowdown |
|---|---|---|---|
| Fortran (DaliLite.v5) | 7.76s | 0.39s | 1.0x |
| Rust (PyO3, release) | 25.37s | 1.27s | 3.3x |

### WOLF-only path (20 pairs)

| Implementation | Time |
|---|---|
| Rust (PyO3, release) | 0.43s |

WOLF is fast (~21ms/pair). The bottleneck is PARSI (branch-and-bound) which
dominates the full pipeline. WOLF+DP+DALICON takes ~0.43s for 20 pairs
vs ~25s total, so PARSI accounts for ~96% of runtime.

## Optimization History

| Change | single_pair | 20_pairs | per_pair | Speedup |
|---|---|---|---|---|
| Baseline | 0.83s | 25.37s | 1.27s | 1.0x |
| LTO+codegen-units=1, HashMap→Vec, stack scoring table | 0.748s | 22.72s | 1.14s | 1.12x |
| SearchState flat ni/ci arrays (eliminate Vec<Vec<i32>> round-trip) | 0.745s | 22.78s | 1.14s | 1.11x |
| target-cpu=native + #[inline(always)] on hot scoring functions | 0.665s | 20.73s | 1.04s | 1.22x |
| + PGO (profile-guided optimization) | 0.646s | 20.0s | 1.00s | 1.27x |
| Bit-packed SearchState (no PGO) | 0.662s | 20.0s | 1.00s | 1.27x |
| **Bit-packed SearchState + PGO** | **0.631s** | **19.0s** | **0.95s** | **1.34x** |

### Round 1: Low-hanging fruit (2026-02-14)

Changes applied:
1. **Compiler**: LTO="fat", codegen-units=1 in workspace release profile
2. **WOLF compare.rs**: HashMap vote counting → flat Vec array (eliminates hashing in inner loop)
3. **PARSI scoring.rs**: `vec![0i32; 101*101]` heap alloc → `[0i32; 101*101]` stack array (hot inner loop)
4. **PARSI stack.rs**: `HashMap<usize, Vec<i32>>` in SearchState → `Vec<Vec<i32>>` (eliminates hashing for all branch-and-bound states)
5. **PARSI stack.rs**: `clearstack` uses `drain()` instead of `.to_vec()` (avoids reallocation)

Result: 12% overall speedup (22.72s vs 25.37s for 20 pairs).
Still 2.9x slower than Fortran (22.72s vs 7.76s).

### Round 2: BinaryHeap + SearchState flat arrays (2026-02-14)

Changes applied:
1. **PARSI stack.rs**: BinaryHeap (max-heap) for O(log n) push/pop instead of sorted Vec O(n) insert
2. **PARSI stack.rs**: SearchState stores flat `ni: Vec<i32>` + `ci: Vec<i32>` directly instead
   of `Vec<Vec<i32>>` candidates. Eliminates all pack/unpack conversion between representations.
3. **PARSI stack.rs**: `copyandput` pushes ni/ci directly (no Vec<Vec<i32>> construction loop)
4. **PARSI stack.rs**: `declump` uses swap-with-last on flat ci array instead of Vec::remove
5. **PARSI stack.rs**: `getnextbest` reads state.ni/ci by reference (zero allocation on pop)

Result: ~0% additional wall-clock improvement (22.78s vs 22.72s, within noise).
The allocation overhead from Vec<Vec<i32>> was not a measurable bottleneck — the inner-loop
scoring computations (get_ess, get_estimate, update_ex) dominate. Code is cleaner and
allocation-free in the branch-and-bound loop, but the real bottleneck is elsewhere.
Still 2.9x slower than Fortran (22.78s vs 7.76s).

### Round 3: target-cpu=native + inline annotations (2026-02-14)

Changes applied:
1. **Compiler**: Added `.cargo/config.toml` with `rustflags = ["-C", "target-cpu=native"]`
   - Enables CPU-specific instructions (AVX2, etc.) on leda's Xeon
2. **PARSI scoring.rs**: `#[inline(always)]` on `get_estimate`, `get_ess`, `segsegscore`, `trimtable`
3. **PARSI stack.rs**: `#[inline(always)]` on `re_estimate`, `getemax_update`

Result: ~9% overall speedup (20.73s vs 22.78s for 20 pairs, consistent across 3 runs).
Single pair: 0.665s (down from 0.745s, 10.7% faster).
WOLF path: 0.320s (down from 0.365s, 12.3% faster).
Now 2.67x slower than Fortran (20.73s vs 7.76s).

### Round 3b: PGO (profile-guided optimization) (2026-02-14)

Steps:
1. Build instrumented: `RUSTFLAGS="-C target-cpu=native -Cprofile-generate=/tmp/pgo-data" maturin develop --release`
2. Run bench.py to collect profile data (24 .profraw files, 6.5MB merged)
3. Merge: `llvm-profdata merge -o merged.profdata *.profraw`
4. Rebuild optimized: `RUSTFLAGS="-C target-cpu=native -Cprofile-use=/tmp/pgo-data/merged.profdata" maturin develop --release`

Result: ~3.5% additional speedup over target-cpu+inline alone (20.0s vs 20.73s for 20 pairs).
Single pair: 0.646s (down from 0.665s). WOLF path: 0.302s (down from 0.320s).
Cumulative from baseline: **21% speedup** (20.0s vs 25.37s).
Now 2.58x slower than Fortran (20.0s vs 7.76s).

Note: PGO requires the profdata file at build time. For reproducible builds, the profdata
should be checked in or regenerated. The `.cargo/config.toml` only has target-cpu=native;
PGO must be applied via RUSTFLAGS manually.

### Round 4: Bit-packed SearchState (2026-02-14)

Changes applied:
1. **PARSI stack.rs**: SearchState stores `bits: Vec<u32>` instead of `ni: Vec<i32>` + `ci: Vec<i32>`.
   Each candidate per segment is 1 bit. State size drops from ~32KB (MAXRES0 * nseg i32s) to
   ~20 bytes (mi[seg] bits per segment, packed into u32 words).
2. **PARSI stack.rs**: PriorityStack precomputes `bit_offset` and `words_per_state` from mi/seglist.
   Added `pack()` and `unpack()` methods for converting between flat ni/ci and bit representation.
3. **PARSI stack.rs**: `copyandput` modifies ni/ci in-place with save/restore (no 32KB clone).
   Only packs the modified state as bits on push.
4. **PARSI stack.rs**: `declump` uses scratch buffers, unpacks/modifies/repacks each state.
5. **PARSI align.rs**: All SearchState construction sites use `pstack.pack()` instead of struct literals.

Result without PGO: 0.662s single, 20.0s 20-pair (~0% over non-bitpack baseline).
Result with PGO: 0.631s single (2.3% over PGO baseline), 19.0s 20-pair (5% over PGO baseline).
Cumulative from baseline: **25% speedup** (19.0s vs 25.37s).
Now 2.45x slower than Fortran (19.0s vs 7.76s).

Analysis: The 2000x state size reduction marginally improves cache behavior under PGO but
does not fundamentally change the performance profile. The remaining 2.45x gap vs Fortran is
dominated by the scoring inner loop (get_ess/get_estimate/segsegscore), where Fortran benefits
from its Fortran-native column-major array layout and aggressive scalar optimization by gfortran.
Further gains would require algorithmic changes (pruning, SIMD scoring) rather than data structure
optimization.
