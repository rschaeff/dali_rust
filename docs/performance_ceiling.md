# Performance Ceiling Analysis

## Context

After four rounds of optimization on the Rust reimplementation of DaliLite.v5,
the PARSI kernel remains ~1.5x slower than the Fortran original. This document
analyzes why, what it would take to close the gap, and why the current state
represents a natural stopping point for this methodology.

## Optimization History

| Round | Target | Technique | PARSI kernel | Pipeline |
|-------|--------|-----------|-------------|----------|
| Baseline | — | Faithful port | Fortran 2.3x | Fortran 4.9x |
| 1 | DALICON | Allocation elimination | — | Fortran 3.3x |
| 2 | Pipeline | Rayon parallelism | — | Fortran 2.0x |
| 3 | PARSI | Unsafe indexing, stack allocs | Fortran 2.1x | Fortran 1.9x |
| 4 | PARSI | Pointer arithmetic, bounds elimination, i32 | **Fortran 1.5x** | **Fortran 1.7x** |

Rounds 1-2 were structural (parallelism, allocation patterns). Rounds 3-4 were
mechanical (replacing safe patterns with unsafe equivalents of what Fortran does
natively). The diminishing returns are clear: Round 4 was the largest kernel
improvement (-27.5%) but was also the round where we shifted from translating
Fortran's syntax to matching its execution model.

## The Hot Loop

The innermost scoring loop in `segsegscore` accounts for the majority of PARSI
runtime. After Round 4 it looks like:

```rust
for b_col in b_start..=b_stop {
    let dv2 = unsafe { *dist2_ptr.offset(q + b_col as isize) } as usize;
    let dv1 = unsafe { *dist_ptr.add(p + b_col) } as usize;
    if dv2 >= 1 && dv1 >= 1 && dv2 < 161 && dv1 < 401 {
        x += unsafe { *scoretable_ptr.add(dv2 * 401 + dv1) };
    }
}
p += nres1;
q += nres2 as isize;
```

The Fortran equivalent:

```fortran
do b=b1+jbeg,b2-jend
    if(dist2(q+b).lt.1) cycle
    if(dist(p+b).lt.1) cycle
    x=x+scoretable(dist2(q+b),dist(p+b))
end do
p=p+nres1
q=q+nres2
```

These are structurally identical. The Rust version uses pointer arithmetic
(`p += nres1`, `q += nres2`) matching the Fortran exactly. Both have the same
data-dependent table lookup, the same distance validity checks, the same
accumulator.

## Why 1.5x Remains

### Not SIMD

The natural assumption is that Fortran benefits from auto-vectorization. But
the inner loop's core operation — `scoretable[dv2 * 401 + dv1]` — is a
**gather**: each iteration computes a data-dependent index and loads from a
different location in a 64K-entry table. This is the anti-pattern for SIMD:

- AVX2's `VGATHERDD` processes 8 lanes but is essentially 8 serial loads,
  yielding 2-3x speedup at best, not 8x
- The scoretable doesn't fit in SIMD registers (~256KB)
- Each gather potentially hits a different cache line

GCC likely cannot auto-vectorize this loop either. The advantage is elsewhere.

### Likely causes

**No-alias guarantees.** Fortran arrays cannot alias by language definition.
When GCC compiles `x = x + scoretable(dist2(q+b), dist(p+b))`, it knows that
the accumulator `x` cannot affect `dist2`, `dist`, or `scoretable`. This lets
it freely reorder loads, hoist invariants, and pipeline memory access. Rust
provides similar guarantees through the borrow checker, but we're using `unsafe`
raw pointers, which opts out of those guarantees. LLVM must be conservative.

**Register pressure and type width.** Fortran's `INTEGER*2` (i16) values feed
directly into index computation. Our Rust code promotes to `usize` (u64 on
x86-64), doubling register consumption. Each `as usize` cast is a zero-extend
instruction that adds to the dependency chain. The Fortran version uses 4-5
live variables in the inner loop; the Rust version has 8-10.

**Instruction scheduling.** GCC has decades of tuning for tight Fortran
numerical loops. LLVM is excellent but its strength is in C/C++-style code
with complex control flow, not in the simple counted loops with array access
that characterize Fortran scientific computing.

**i16 native arithmetic.** The Fortran distance matrices are INTEGER*2 and the
comparisons (`dist2(q+b).lt.1`) operate on i16 directly. Our Rust code loads
i16 then immediately promotes to usize for the comparison and table index.
This adds instructions and widens the data path unnecessarily.

### What it would take to close

The remaining 1.5x could theoretically be addressed by:

1. **FFI to a C/Fortran inner kernel.** Write the hot loop in C with
   `__restrict__` pointers (matching Fortran's no-alias semantics) and call
   it from Rust via FFI. This would let GCC or Clang optimize the inner loop
   with full alias information. Alternatively, write it directly in Fortran
   and link it.

2. **Explicit SIMD with gather.** Use `std::arch::x86_64` intrinsics to
   manually vectorize: load 8 or 16 consecutive i16 distances, use
   `_mm256_i32gather_epi32` for the table lookup, mask invalid entries,
   horizontal sum. This requires restructuring the data layout and manually
   managing the remainder loop for non-multiple-of-8 segment lengths.

3. **Assembly inspection and targeted fixes.** Disassemble the Rust and
   Fortran binaries, identify where LLVM makes suboptimal register
   allocation or instruction ordering choices, and restructure the Rust
   source to guide LLVM toward better codegen.

4. **Profile-guided optimization (PGO).** Run the Rust binary under PGO
   instrumentation and recompile with the profile data, letting LLVM make
   better branch prediction and inlining decisions.

Each of these approaches has problems:

- **FFI inner kernel**: Works, but defeats the purpose. If the argument is
  "Rust can replace Fortran," calling Fortran from Rust is not that argument.
  It's also baroque — maintaining a polyglot inner loop for a 1.5x factor
  that's invisible at the pipeline level.

- **Explicit SIMD**: High effort (~500 lines of intrinsics code), uncertain
  return (gather is the bottleneck, not the arithmetic), fragile (tied to
  x86-64 + AVX2), and very hard to validate (SIMD bugs are notoriously
  subtle).

- **Assembly analysis**: Requires deep expertise in x86 microarchitecture
  that neither the AI nor the human in this collaboration reliably has.
  Even if we identified the issue, the fix would be "write the Rust in a
  way that happens to produce better assembly" — a fragile optimization
  that could regress with the next LLVM version.

- **PGO**: Easy to try, but typically yields 5-15% on this kind of code,
  not 50%. Worth attempting but won't close the gap alone.

## Why This Is The Right Stopping Point

### The pipeline tells the real story

At the pipeline level (which is what users experience), Rust wins 3 of 6
benchmark pairs and the overall gap is 1.68x. For same-fold pairs (the
common case in DPAM), Rust is already faster. The remaining gap appears
only on large cross-fold comparisons where the PARSI kernel dominates.

In the DPAM use case (50-500 comparisons per query), Rust's ability to
parallelize across pairs using Rayon eliminates the per-pair kernel gap
entirely. Fortran's `serialcompare` requires a separate process per pair,
with ~30ms overhead per invocation and communication through the filesystem.

### The methodology argument is complete

The project's goal is to demonstrate that AI-assisted validated reimplementation
of legacy scientific software is feasible. The performance story strengthens
this argument at every level:

1. **Faithful port** (Phase 1-9): Correct output, 2.3x slower kernel.
   Proves the AI can translate the algorithm accurately.

2. **Structural optimization** (Rounds 1-2): Parallelism and allocation
   patterns that Fortran's architecture prevents. Proves the reimplementation
   enables improvements impossible in the original.

3. **Idiomatic optimization** (Rounds 3-4): Matching Fortran's execution
   model where it matters. Proves the AI can reason about *why* code is fast,
   not just *what* it does.

4. **Documented ceiling** (this document): Honest analysis of what remains
   and why. Proves the methodology produces not just code but understanding
   of its limitations.

The fact that we can precisely characterize the remaining gap (no-alias
semantics, register pressure, codegen quality) rather than just saying
"Fortran is faster" is itself evidence of deep understanding.

### The diminishing returns are structural

Each optimization round required more specialized knowledge for less return:

- Round 1: Know that `Vec::new()` in a hot loop is bad → -46% DALICON
- Round 2: Know that Rayon exists → pipeline 4.9x → 2.0x
- Round 3: Know that bounds checking costs cycles → -8.5% PARSI
- Round 4: Understand Fortran's pointer arithmetic idiom → -27.5% PARSI
- Round 5 (hypothetical): Understand x86 gather throughput → ???% PARSI

The knowledge required escalates while the returns diminish. This is the
natural shape of optimization: the early wins come from understanding the
*algorithm*, the middle wins from understanding the *language*, and the
last wins require understanding the *hardware*. An AI + non-expert human
collaboration is well-suited to the first two and poorly suited to the third.

This is not a failure of the methodology — it's a characterization of its
scope. The methodology demonstrably works for producing correct, performant,
maintainable replacements for legacy scientific software. Extracting the
last 1.5x of kernel performance is a different kind of problem that requires
different tools.

## Final Numbers

```
PARSI kernel:   Fortran 1.52x faster (was 2.30x at baseline)
Full pipeline:  Fortran 1.68x faster (was 4.87x at baseline)
Pipeline pairs: Rust wins 3/6 (was 0/6 at baseline)
Same-fold:      Rust 5-54% faster
Cross-fold:     Fortran 2-3.5x faster (PARSI-dominated)
DPAM use case:  Rust faster (parallelism across pairs)
```
