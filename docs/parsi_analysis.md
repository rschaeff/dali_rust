# PARSI Module Deep Analysis

## Overview

PARSI is an exhaustive branch-and-bound alignment optimizer. It takes the DALICON output (rough alignments) and systematically searches for the optimal structural alignment by:
1. Decomposing the protein into a **domain tree** of SSE segments
2. Scoring all possible segment-to-residue mappings using a **DaliLite scoring function**
3. Using **branch-and-bound** with a priority queue (scored stack) to enumerate alignments from best to worst
4. **Splitting** the search space at each step by dividing candidate residues for the most informative segment
5. Outputting multiple sub-optimal alignments per domain

## Source Files (4765 lines total)

| File | Lines | Purpose |
|------|-------|---------|
| `parsi-1.f` | 311 | Entry-level routines: `perresiduescore`, `findhomolog`, `lsave`, `singletstack`, `ltransversion` |
| `parsi-admin.f` | 1095 | Administration: protein I/O, distance matrices, scoring tables, domain tree, search space init |
| `parsi-align.f` | 1330 | Core alignment: `align()` orchestrator, `dostack`, `loadstack`, `compress`, `shipup`, `saveit`, `output1/2` |
| `parsi-score.f` | 811 | Score calculation: `segsegscore`, `singletex`, `doubletex`, `update_ex`, `get_ess`, `get_estimate`, `trimtable` |
| `parsi-stack.f` | 1218 | Stack/priority queue: `split`, `getnextbest`, `putrequest`, `toprequest`, `declump`, `clearstack`, bit-packing, linked list management |

## Constants (from parsizes.for)

```fortran
maxres1 = 4000      ! max residues in protein 1 (compressed)
maxres2 = 4000      ! max residues in protein 2
maxres0 = 400       ! maxres2/10 — max candidates per segment (10-block stepping)
exdim = 1000000     ! size of score table ex()
maxstack = 10000    ! stack size (in integers)
maxstack2 = 10000   ! link array size
maxseg = 200        ! max SSE segments
maxdom = 400        ! max domain tree nodes (maxseg*2)
infinit = 10000000  ! sentinel for "infinity"
nul = -99           ! sentinel for "unaligned"
boxdim = 501        ! link_pointer grid dimension for fast insertion
```

## Key Global State

The `parsi` module (in comparemodules.f) maintains persistent state between calls:
- `dist(maxres1*maxres1)` — int16 distance matrix for protein 1 (compressed CA coordinates)
- `dist2(maxres2*maxres2)` — int16 distance matrix for protein 2
- `dist1sum(maxres1,0:maxres1)` — cumulative distance sums for protein 1
- `dist2sum(maxres2,0:maxres2)` — cumulative distance sums for protein 2
- `segmentrange(2,maxseg)` — compressed SSE ranges for protein 1
- `segmentrange0(2,maxseg)` — original SSE ranges for protein 1
- `ss(maxseg,maxseg)` — self-score matrix (segment-segment scores)
- `upper/lower(maxseg,maxseg)` — distance sum bounds (70%/130% of self-distances)
- `ndom, node_child(2,maxdom), domns(maxdom), domseglist(maxseg,maxdom)` — domain tree
- `lfix(maxseg,maxdom)` — which segments are "fixed" in each domain node
- `cut(maxdom)` — score cutoffs per domain node
- `weight(1000)` — Gaussian envelope weights
- `scoretable(160,400)` — precomputed scoring table

Module parameters:
- `lpreali = .false.` — no prealignment exclusion
- `lnul = .true.` — allow unaligned segments
- `ldb = .true.` — read from list (batch mode)
- `lseqtl = .true.` — enforce sequential topology

## Data Flow

### Initialization (`init_parsi`)
1. `weights()` — Gaussian envelope: `weight[i] = nint(100*exp(-(i/20)^2))`, clipped to 0 if <5
2. `fillscoretable()` — 160×400 table: `scoretable[b,a] = nint(100*weight[a]*(0.20-abs(x-y)/x))`
   - where `a` indexes distance in protein 1 (a/10.0 Å), `b` indexes binned distance in protein 2
   - Binning: b<100 → y=b*0.1; 100≤b<125 → y=10+(b-100)*0.4; b≥125 → y=b-125+20

### Per-pair: `dowork_parsi(cd1, cd2, outputunit, oldcd1, datpath1, datpath2, lfirstonly)`

#### Phase 1: Load protein 1 (cached if cd1 unchanged)
`getstructure1()`:
1. `getprotein1()` → `parsireadproteindata()`: Read .dat file (same format as WOLF/DP)
2. `treehack()`: Split multi-segment leaf nodes into binary tree
3. `compressca()`: Remap CA coords to packed indices (remove inter-segment gaps)
4. `getdist()`: CA distance matrix for protein 1, stored as int16 (distance × 10, capped at 1000)
5. `selfscore()`: Segment-segment self-scores using `weight[]` (Gaussian envelope × 20 per residue pair)
6. `getupperlower()`: 70%/130% distance-sum bounds per segment pair
7. `hackdist()`: Cap distances at 400 (= 40.0 Å)
8. `flex()`: Determine fixed/free segments per domain node (density-based, cutoffs 0.70 and 0.90)
9. `setcut()`: Score cutoffs per domain: `9000*(polynomial in lali)`, capped at 200 residues
10. `getdist2sum()`: Cumulative distance sums for protein 1 (int4: max(100, nint(1000*dist)))
11. `setldom()`: Mark domains with ≥`startsize` segments as "align"
12. `setminseglen()`: Min segment length: 6 for E, 8 for H
13. `setngap()`: N-terminal gap positions per segment

#### Phase 2: Load protein 2
`getsequence2()`:
1. `getprotein2()` → `parsireadproteindata()`: Read .dat file
2. Build `segment2[]`: residue-to-segment mapping
3. `getdist2sum()`: Cumulative distance sums for protein 2
4. `getdist2()`: Binned distance matrix for protein 2:
   - dist ≤ 10 Å: `nint(dist/0.1)` (bins 1-100)
   - 10 < dist ≤ 20: `100 + nint((dist-10)/0.4)` (bins 100-125)
   - 20 < dist < 55: `125 + nint(dist-20)` (bins 125-160)
   - dist ≥ 55: 160 (default)

#### Phase 3: Alignment (`align()`)

This is the core algorithm. It processes the domain tree bottom-up.

##### Search space representation:
- `trans(ir, iseg)` — maps candidate index `ir` for segment `iseg` to a residue number in protein 2 (or `nul`)
- `ci(ir, iseg)` — candidate indices (subset of 1..mi)
- `ni(iseg)` — number of active candidates for segment `iseg`
- `mi(iseg)` — maximum candidates (full search space)
- `ex(exdim)` — flat array of precomputed scores:
  - `ex[start(iseg,iseg) + ir]` = singlet score for segment `iseg` at candidate `ir`
  - `ex[start(aseg,bseg) + (ir-1)*mi(bseg) + jr]` = doublet score for segments `aseg`,`bseg` at candidates `ir`,`jr` (where aseg<bseg)

##### `init_searchspace()`:
- For each segment, enumerate candidate residues in protein 2 at 10-residue steps (bl=10)
- Add N-terminal gap residues (negative indices down to -29)
- Add `nul` candidate (unaligned)
- If `lpreali=true`: read `.cutz` prealignment file and exclude over-aligned positions

##### `update_ex()` — Score calculation:
For each segment pair not yet scored:
1. **Singlets** (`singletex`): For each candidate residue for segment `iseg`:
   - Try all 10 shift positions (bl=10) within the 10-block
   - For each shift, compute per-residue-pair score using `scoretable[dist2[q], dist[p]]`
   - Trim destabilizing end residues (`trimtable`) — iteratively remove N/C-terminal residues with negative row sums until min segment length reached
   - Remember best score across shifts
   - Store `ibeg/iend` trimming offsets in `s_beg/s_end` arrays

2. **Doublets** (`doubletex`): For each pair of candidates (ir, jr) for segments (aseg, bseg):
   - Call `segsegscore()`: tries all ishift×jshift combinations (0..bl-1 each)
   - Distance sum filter: quick reject if cumulative distance sum outside 70-130% of self
   - Precise calculation: sum `scoretable[dist2[q+b], dist[p+b]]` over aligned residue pairs
   - Enforce sequential topology (`lseqtl`): if iseg < jseg, require transires < transjres
   - Disallow segment overlaps in protein 2

##### `align()` main loop (bottom-up over domain tree):
```
for idom = ndom down to 1:
    if not ldom(idom): skip

    shipup(): build work search space from child node results
    compress(): remap ex[] indices after merging child nodes
    update_ex(): compute scores for new segment pairs

    if domns(idom) == 1:
        output1(): enumerate all singlet candidates above cutoff
    elif domns(idom) == 2:
        output2(): enumerate all doublet combinations above cutoff
    else:
        loadstack(): initialize priority queue from child alignments

        loop:
            dostack() → getnextbest():
                pop best from priority queue
                if unique (all ni==1): return alignment (flag=5)
                if empty/below cutoff: return (flag=1/2)
                if overflow: clearstack (flag=3)
                else: split() search space and push halves

            if flag != 5: break
            declump(): remove found alignment from stack
            if degenerate: continue (skip duplicates)
            findhomolog(): find matching segments in child alignments
            perresiduescore(): refine alignment at single-residue resolution
            saveit(): store alignment for parent nodes
            segmentpairs(): output segment pair summary

            if lfirstonly and icyc >= 5: break
```

##### `split()` — Search space division:
1. Find segment `xseg` with highest score contribution (via `get_ess`)
2. Compute `emax[xres]` = max score achievable from each candidate for xseg
3. Split candidates into top 10% (ci1) and bottom 90% (ci2) by emax
4. If sequential topology: further split top half into left/right/nul groups
5. Push each subgroup onto priority queue with re-estimated scores

##### `perresiduescore()` — Single-residue refinement:
- Takes the 10-block alignment and refines to individual residue positions
- For each segment, enumerate residues in a ±9 window around the 10-block position
- Run full branch-and-bound at residue resolution (reusing the same stack infrastructure)

##### Stack/Priority Queue:
- Stack entries are bit-packed: `chunk = 1 + ceil(sum(mi)/31)` integers
  - First integer = score estimate
  - Remaining integers = bit-packed candidate masks (1 bit per candidate per segment)
- Linked list maintains score ordering: `link(1,i)` = next lower, `link(2,i)` = next higher
- `link_pointer[box]` = fast lookup grid indexed by `int(est/firstest * boxdim)`
- `bestest` always points to the highest-scoring entry

### Output Format

```
refine<cd1><cd2> <idom> <score> <nseg> <a1_1> <a2_1> ... <a1_n> <a2_n> <b1_1> <b2_1> ... <b1_n> <b2_n>
```

Where:
- `cd1cd2` = concatenated 5-char structure codes
- `idom` = domain tree node index (can be 1, 2, or 5 for different tree levels)
- `score` = alignment score (from the scoring function)
- `nseg` = number of aligned segment pairs
- `a1_i, a2_i` = start, end of segment i in protein 1 (original residue numbering via `segmentrange0`)
- `b1_i, b2_i` = start, end of segment i in protein 2
- Unaligned segments: `a1 = a2 = b1 = b2 = nul = -99`

## Complexity Analysis

For a protein pair with `nseg` segments and `nres2` residues in protein 2:
- Candidates per segment: ~nres2/10 + ngap + 1(nul) ≈ nres2/10
- Score table `ex[]`: ~O(nseg² × (nres2/10)²) entries
- Branch-and-bound: exponential worst case, but pruning makes it practical
- Per-residue refinement adds another factor but operates on small neighborhoods

## Key Differences from WOLF/DP/DALICON

1. **Domain tree**: PARSI processes a hierarchical decomposition, not flat pairwise
2. **10-block stepping**: Candidates are at 10-residue intervals, then refined to single residues
3. **Bit-packed stack**: Search states are compressed using bit packing to save memory
4. **Two-phase search**: Coarse (10-block) then fine (single residue) resolution
5. **Multiple outputs**: Generates multiple sub-optimal alignments per domain, not just the best
6. **Sequential topology enforcement**: `lseqtl=true` requires N-to-C ordering
7. **Compressed coordinates**: Protein 1 CA coords are repacked to eliminate inter-segment gaps

## Shared Code with Existing Modules

PARSI does NOT reuse WOLF/DP/DALICON code. It has its own:
- Distance matrix computation (`getdist`, `getdist2` — different binning from DP)
- Scoring function (`segsegscore` — Dali-like but with different parameters)
- Protein reader (`parsireadproteindata` — reads all 4 sections of .dat including domain tree)
- Superposition is NOT used (PARSI operates purely in distance space)

However, the **conceptual** scoring (Gaussian-weighted distance comparison) is similar to DP's scoring.

## Implementation Considerations

### Most Complex Parts
1. **Stack/priority queue with bit-packing** (~500 lines): The linked list with box-indexed fast lookup, bit-packed candidate masks, and swap/remove operations is the most intricate part
2. **align() orchestrator** (~330 lines): Bottom-up domain tree traversal with compress/shipup/loadstack
3. **segsegscore** (~100 lines): The inner scoring loop with 10-block shifts, trimming, distance sum filters, and sequential topology checks

### Moderate Complexity
4. **split()** (~240 lines): Dividing search space by score percentile
5. **update_ex/singletex/doubletex** (~230 lines): Populating the score table
6. **perresiduescore** (~230 lines): Single-residue refinement

### Straightforward
7. **Protein I/O and setup** (~400 lines): Reading .dat files, building distance matrices, domain tree
8. **Utility functions** (~200 lines): setcut, flex, treehack, etc.

### Potential Pitfalls
1. **Fortran 1-indexing in bit-packing**: `bitpack`/`bitunpack` use bits 1-31, `ishft` operations
2. **Compressed coordinate mapping**: `compressca` repacks residues, `segmentrange` vs `segmentrange0`
3. **ex[] flat array layout**: Must match Fortran's `start(iseg,jseg)` + offset scheme exactly
4. **Distance matrix format difference**: `getdist` (int16, dist×10, cap 1000) vs `getdist2` (binned: 0.1/0.4/1.0 Å bins)
5. **Domain tree bottom-up traversal with state accumulation**: `ali_save`, `ali_start`, `compress`, `shipup`
6. **`getemax` uses `ran(seed)`**: Another instance of the gfortran intrinsic RAND shadowing issue (seed=1234567, effectively constant)
7. **Stack overflow handling**: `clearstack` removes bottom 90% on overflow
8. **`treehack` modifies domain tree**: Splits multi-segment leaf nodes, increasing ndom

## Estimated Implementation Size

~2000-3000 lines Python, organized as:
```
parsi/
├── __init__.py          # Public API
├── parsizes.py          # Constants
├── protein_io.py        # parsireadproteindata, compressca, getdist, getdist2
├── scoring.py           # weights, fillscoretable, segsegscore, singletex, doubletex
├── domain_tree.py       # treehack, flex, setcut, setldom, domain tree traversal
├── search_space.py      # init_searchspace, update_ex, get_ess, get_estimate
├── stack.py             # Priority queue: putrequest, toprequest, getnextbest, split, declump
├── bitpack.py           # bitpack, bitunpack, putchunk, getchunk, translatechunk
├── align.py             # align(), dostack, loadstack, compress, shipup, perresiduescore
└── parsi.py             # dowork_parsi orchestration
```
