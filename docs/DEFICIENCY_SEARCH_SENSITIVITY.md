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

## Where to look

Two observations that should narrow it, from this data rather than from inspection:

1. **The misses start far above both cutoffs.** `ZCUT_DEFAULT = 2.0` (`dali-core/src/dp/`)
   and `zcut1 = 1.0` (`dali-core/src/filter95/`) cannot explain missing 19% of templates
   whose Fortran Z is >= 8. A cutoff *may* account for part of the Z 2-6 band, where a
   genuine Rust alignment could score under 2.0 and be discarded. **First experiment: rerun
   a protein set with both cutoffs at 0 and re-measure recall by band.** That cleanly splits
   thresholding from search, and it is cheap.

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
