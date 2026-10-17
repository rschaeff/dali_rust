# Deficiency Report: Search Sensitivity at ECOD70 Scale

**Date:** 2026-10-13
**Tested against:** DPAM v2.0 step 07 (Iterative DALI) through step 24 (Integrate Results)
**Test set:** 9,087 *Rhizomucor* proteins (full Stage B baseline), 48,385 distinct ECOD70
templates, 113,120 Fortran alignments
**Reference backend:** Fortran DaliLite.v5, via DPAM's own `run_dali` per-template path
**Harnesses:** `~/work/fungi/analysis/dali_bench/` (bench_dali_v3.py, e2e/)
**Severity:** BLOCKING for the DPAM use case
**Status:** OPEN

## Executive Summary

`DEFICIENCY_REPORT_DPAM_INTEGRATION.md` (2026-02-15) closed out four integration
deficiencies on two proteins and attributed the residual hit-count gap to "inherent search
divergence in WOLF/PARSI heuristics (documented at 78-87% match rate)." **At ECOD70 scale
that gap is not 13-22%, it is 63%, and it changes DPAM's final answer.**

| Metric | n=2 (Feb report) | n=9,087 (this report) |
|--------|------------------|----------------------|
| Template recall | 0.82 / 0.49 | **0.375** |
| Recall at Fortran Z>=8 | not measured | 0.808 |
| Rust-only hits | 0 | 803 domains (vs 7,245 lost) |
| Alignments | 124 of 151, 28 of 57 | 32,471 of 113,120 |

End-to-end, with steps 1-6 restaged from the Fortran baseline, step 7 swapped, and steps
08-24 re-run under the Fortran backend in **both** arms so step 7 is the only variable:

    base domains retained positionally   7,981 / 15,226  (0.524)
    proteins losing EVERY domain         1,547 / 9,087   (0.170)
    proteins with a different count      4,105 / 9,087   (0.452)
    residue Jaccard, mean                0.688
    of domains that survive, T-group     7,308 / 7,981   (0.916)

**One protein in six loses its entire domain call.** The 92% T-group agreement that looks
reassuring is conditional on surviving at all — it describes the surviving half.

The validity gate passed: restaging with Fortran's own step-7 hits reproduces the baseline
domain tables **300/300 byte-identical**, so none of the above is staging loss.

## The deficit amplifies downstream

This is the part that makes it blocking rather than merely a quality gap. Take the ECOD
domain each baseline call selected, look up Fortran's best Z for that template, and split by
whether the Rust arm reproduced the call:

    band      kept    lost   retained
    Z>=8     6,249   3,572     0.636
    6<=Z<8     597   1,250     0.323
    4<=Z<6     351   1,870     0.158
    Z<4        111   1,226     0.083

Template recall at Z>=8 is 0.808, but only **0.636 of domain calls resting on a Z>=8
template survive**. DPAM's domain judge reads `Hcount`, `Scount` and `Hit_cov` — counts of
*supporting* alignments — so dropping two thirds of the alignments changes the verdict even
where the top template was recovered. A per-template recall figure therefore understates
the damage; the twilight loss does not stay in the twilight.

## What is excluded as the cause

Each of these was measured, not argued:

- **DSSP divergence (Deficiency 2) is not it, and the `.dat` bypass does not fix this.**
  Feeding Rust a Fortran-provenance `.dat` cache built with `import.pl` makes recall
  slightly *worse*: 0.307 vs 0.329 paired, 0.922 vs 0.955 at Z>=8, while the rate holds
  (15.55x paired). The mitigation in the Feb report addresses the representation, and the
  representation was not the limiter. The SSE divergence is real — 78% of 198 sampled
  templates differ, usually an off-by-one boundary — it just is not what loses the hits.
- **`.dat` cache coverage.** A no-cache run and a rust-`.dat` run produce byte-identical
  hit lists, so `dat_dir` is purely a speed cache with no effect on sensitivity.
- **Template import failures.** All 63,065 ECOD70 templates were run through
  `dali.import_pdb`: 172 fail (0.27%), all `TooFewResidues`. The `TooFewResidues(0)` stream
  in the logs is loud and immaterial.

By elimination the deficit is in the alignment search itself.

## Measured: the cutoffs are not the problem (2026-10-13, jobs 685601/685602)

The first experiment proposed below has been run. Both cutoffs to 0 over 300 proteins and
190,482 template calls, against the same `dat_cache`, pool width and candidate lists, with
a default-cutoff control arm that reproduced the compared-against arm **20/20
byte-identical** so the result could not be the rebuild (`.cargo/config.toml` sets
`target-cpu=native` and the nodes are Xeon Phi 7250).

Every Fortran template falls in exactly one of three classes:

    recalled at default cutoffs    29,333  (0.375)
    THRESHOLDED (zcut0 recovers)    7,334  (0.094)
    NEVER FOUND (search miss)      41,561  (0.531)
    recall: 0.375 -> 0.469

    band      fortran  recalled  thresh   never  rec@def  rec@z0
    Z>=8        22490     18170     417    3903    0.808   0.826
    6<=Z<8       8555      4892     736    2927    0.572   0.658
    4<=Z<6      16870      4663    2838    9369    0.276   0.445
    Z<4         30313      1608    3343   25362    0.053   0.163

**At Z>=8 the cutoffs are worth 1.8 points, 0.808 -> 0.826.** 3,903 confident templates are
never found at all. What the cutoffs do hold back is the twilight — Z 4-6 from 0.276 to
0.445, Z<4 from 0.053 to 0.163 — which is the band DPAM needs least. Opening them is also
not free: templates hit by Rust but not Fortran go from 166 to 927, a 5.6x precision cost,
and recall still only reaches 0.469. **Lowering the cutoffs is not a fix.**

The alignments being discarded are genuine, not noise: median Rust Z 2.63 against Fortran's
4.10, and 81% within 2x of Fortran's Z.

### The gate is the pre-DALICON DP pass, on unrefined scores

All 7,334 recovered templates end with Z >= 2.0, so the FINAL DP gate did not drop any of
them. A second arm with `DALI_ZCUT_DP=0` and FILTER95 left at its default 1.0 then returned
hits **300/300 byte-identical to the both-cutoffs arm**, which locates it exactly:

- **FILTER95/pipe96 at zcut 1.0 discards nothing.** Not a contributor. (`PIPENBEST = 2` may
  still bind on the PARSI path; `DALI_PIPENBEST` is there to test that separately.)
- **The whole thresholded population dies at the first `run_dp` call**, before DALICON.
  That pass scores *unrefined* WOLF output at the same 2.0 as the post-refinement pass.
  DALICON refinement raises scores — these hits end up at a median Z of 2.63 once refined —
  so a single threshold applied before refinement discards alignments that would clear it
  after. Fortran orders its pipeline the same way, so this is not by itself a divergence
  from Fortran; it is where our losses are concentrated, and the cheapest knob available.
  Whether Fortran's pre-DALICON scores sit higher than ours on the same candidates is the
  obvious next question, and would separate "same gate, worse scores" from "same scores,
  fewer candidates".

Net: thresholding is real but second-order and twilight-only. **Candidate generation carries
53% of the deficit, including 3,903 templates at Z>=8, and that is the deficiency.** The
second pointer below is now the live one.

## Where to look

Two observations that should narrow it, from this data rather than from inspection:

1. **The misses start far above both cutoffs.** `ZCUT_DEFAULT = 2.0` (`dali-core/src/dp/`)
   and `zcut1 = 1.0` (`dali-core/src/filter95/`) cannot explain missing 19% of templates
   whose Fortran Z is >= 8. A cutoff *may* account for part of the Z 2-6 band, where a
   genuine Rust alignment could score under 2.0 and be discarded. This was measured and is answered above: the cutoffs account for
   9.4 points of recall, 1.8 of them at Z>=8.

2. **The asymmetry points upstream of PARSI refinement.** `parsi_mismatch_analysis.md`
   found Python/Rust sometimes finding *genuinely higher-scoring* alignments than Fortran
   (1a87A vs 1a00A: 2,227,558 vs 1,140,652, independently verified). Yet at scale Rust's
   hits are a near-**subset** of Fortran's — 803 Rust-only domains against 7,245 lost. A
   refinement that is sometimes better but a hit list that is almost strictly smaller is not
   a refinement-quality problem; it is candidate generation. That points at WOLF (spatial
   hash, `NEIBORCUTOFF`, segment-pair enumeration) and at how many DCCP candidates survive
   to refinement, not at PARSI's branch-and-bound quality.

The "78-87% match rate" in module validation is consistent with all of this and does not
contradict it: those tests measure agreement on pairs that both backends align. They cannot
see a template that never becomes a candidate, which is exactly the failure mode here.

## Rate, for the record

The speed is real and was measured paired, on the 398 templates where both backends emitted
an alignment under identical 272-way pool contention: **14.3x per alignment** (paired
cpu-time 8,857.5s vs 592.3s = 14.95x; Fortran emits 1.05x more alignments per shared
template). Pool wall-clock is 7.0x at 406 of 1,367 alignments.

Quote the 14.3x only with "at 0.375 template recall" attached. In DPAM terms the current
trade is 14x throughput for one protein in six losing its domain call, and the decision
there is already made: **DPAM stays on DaliLite.v5.** Closing this deficiency is what would
reopen it.

## Two unrelated blockers in the DPAM integration

Independent of sensitivity, the fast path is not reachable in production today:

- `run_step7` forces `use_batch=True` whenever `DPAM_DALI_BACKEND` is `rust` or `auto`, so
  the per-template `run_dali` path that was benchmarked here is never taken by the pipeline.
  `RustDALI.batch_search` answers a different question — "the globally best match per
  masking round" — and produced 8 alignments against Fortran's 1,367 on the same 6 proteins.
  It must not be benchmarked against step 07.
- `template_dat_dir` has no caller and no CLI flag in `dpam_c2`, so a production Rust run
  would regenerate `.dat` from PDB per protein rather than reuse a shared ECOD70 cache.

## Reproduction

    # paired per-template rate + recall by Z band
    cd ~/work/fungi/analysis/dali_bench
    sbatch run_bench.sh                      # jobs 684842, 685398

    # end-to-end: stage arms, swap step 7, re-run steps 08-24, compare step-24 tables
    sbatch e2e/run_e2e.sh                    # job 685466, n=300 seeded sample
    sbatch --array=1-16 e2e/run_e2e_array.sh # job 685467, full 9,087
    python3 e2e/compare_e2e.py e2e/full.txt

    # importer survey over all of ECOD70
    sbatch e2e/run_probe_imports.sh          # job 685478

Full write-up and the per-domain disagreement table (12,859 rows) are in
`~/work/fungi/analysis/dali_bench/e2e/RESULTS.md` and `e2e/disagreements.tsv`.

Two traps worth knowing if this is re-measured: step-7 hits are keyed by ECOD **uid**
(`000153702`) while step 24's `ECOD_num` is the **domain id** (`e4bpmA1`), so joining them
needs the `ECOD_length` translation or it silently returns zero rows; and
`step24_integrate/` is created only when a table is written, so a protein that legitimately
integrates zero domains leaves no directory and is indistinguishable on the file tree from
a crash.
