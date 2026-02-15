//! DALICON module: Monte Carlo genetic algorithm refinement.
//!
//! Translated from comparemodules.f dowork_dalicon(), lean.f, clean.f, testi.f.
//!
//! Key Fortran compatibility notes:
//! - All scoring arithmetic uses f32 to match Fortran `real` precision
//! - RNG is a constant (gfortran RAND intrinsic quirk): always returns 0.0878167152f32
//! - ca1 coordinates persist across comparisons with same cd1 (DaliconCd1State)
//! - getblocks uses two separate `if` statements, not `else if`
//! - Distance matrices are 1-based (DaliconDistMat)

use ndarray::Array2;

use crate::numerics::fitz::fitz;
use crate::numerics::kabsch::{u3b, transrotate};
use crate::numerics::rounding::nint;
use crate::store::ProteinStore;
use crate::types::{AlignmentBlock, DaliconDistMat};

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

const REFITOL: f32 = 2000.0;
const MAXITER: usize = 100;
const FITZRCUT: f64 = 4.0;
const FITZMAXITER: usize = 10;
const BLOCKSIZE: usize = 4;
const IMAX0: usize = 50;
const ITRIM2: usize = 5;
const MAXPAIR: usize = 6000;
const MAXRES: usize = 5000;

/// gfortran RAND(124462287) constant — see dalicon_notes.md
const RNG_CONSTANT: f32 = 0.0878167152;

// ---------------------------------------------------------------------------
// Scoring (f32 arithmetic to match Fortran)
// ---------------------------------------------------------------------------

/// Metropolis acceptance probability.
#[inline]
fn getp(x: f32) -> f32 {
    if x > 0.0f32 {
        1.0f32
    } else if x < -5.0f32 {
        0.0f32
    } else {
        x.exp()
    }
}

/// Compute Gaussian weight table (f32): w[i] = exp(-0.0025 * i^2).
/// Thread-safe via OnceLock.
fn gagaweights() -> &'static [f32; 101] {
    use std::sync::OnceLock;
    static WEIGHTS: OnceLock<[f32; 101]> = OnceLock::new();
    WEIGHTS.get_or_init(|| {
        let mut w = [0.0f32; 101];
        for i in 0..101 {
            let fi = i as f32;
            w[i] = (-0.0025f32 * fi * fi).exp();
        }
        w
    })
}

/// DALICON scoring function using precomputed weights.
#[inline]
fn scorefun_w(a: i16, b: i16, wght: &[f32; 101]) -> f32 {
    let ia = a as i32;
    let ib = b as i32;
    let x = (ia - ib).unsigned_abs() as f32 / 10.0f32;
    let y = (ia + ib) as f32 / 20.0f32;
    let yint = nint(y as f64);
    if y > 100.0f32 {
        return 0.0f32;
    }
    let iy = yint.clamp(0, 100) as usize;
    let s = if y > 0.0f32 {
        wght[iy] * (0.20f32 - x / y)
    } else {
        wght[iy] * 0.20f32
    };
    s * 100.0f32
}

// ---------------------------------------------------------------------------
// Distance matrix (1-based, for DALICON)
// ---------------------------------------------------------------------------

/// Compute 1-based distance matrix for DALICON.
///
/// d[i,j] = nint(10 * dist(ca[:,i-1], ca[:,j-1])) for i,j in 1..=nres.
/// d[0,:] and d[:,0] are zero padding.
pub fn getgagadist(ca: &Array2<f64>, nres: usize) -> DaliconDistMat {
    let size = nres + 1;
    let mut data = Array2::<i16>::zeros((size, size));
    for i in 1..=nres {
        for j in (i + 1)..=nres {
            let dx = ca[[0, i - 1]] - ca[[0, j - 1]];
            let dy = ca[[1, i - 1]] - ca[[1, j - 1]];
            let dz = ca[[2, i - 1]] - ca[[2, j - 1]];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            let d = nint(10.0 * dist) as i16;
            data[[i, j]] = d;
            data[[j, i]] = d;
        }
    }
    DaliconDistMat { data, nres }
}

// ---------------------------------------------------------------------------
// Testi: tetrapeptide seeding
// ---------------------------------------------------------------------------

/// Compute left/right bounds from prealignment.
/// Translated from ssap.f getleftrite().
fn getleftrite(
    nres1: usize,
    nres2: usize,
    preali1: &[i16],
    width: i32,
) -> (Vec<i32>, Vec<i32>) {
    let size = nres1 + 1;
    let mut left = vec![0i32; size];
    let mut rite = vec![nres2 as i32; size];

    // Forward pass
    for i in 2..size {
        if preali1[i] > 0 {
            left[i] = left[i - 1].max(preali1[i] as i32 - width);
        } else {
            left[i] = left[i - 1];
        }
    }

    // Backward pass
    for i in (1..nres1).rev() {
        if preali1[i] > 0 {
            rite[i] = (preali1[i] as i32 + width).min(rite[i + 1]);
        } else {
            rite[i] = rite[i + 1];
        }
    }

    (left, rite)
}

/// Build tetrapeptide candidate pool from prealignment and distance maps.
///
/// Translated from testi.f testi().
fn testi(
    nres1: usize,
    nres2: usize,
    preali1: &[i16],
    d1: &Array2<i16>,
    d2: &Array2<i16>,
    wght: &[f32; 101],
) -> Vec<(usize, usize)> {
    // Collect aligned positions
    let aligned: Vec<usize> = (1..=nres1)
        .filter(|&i| preali1[i] != 0)
        .collect();

    let (left, rite) = getleftrite(nres1, nres2, preali1, 10);

    // Build score map using HashMap for sparse storage
    let mut score_map = std::collections::HashMap::new();
    for i in 1..=nres1 {
        for j in (left[i] + 1)..rite[i] {
            let j = j as usize;
            if j < 1 || j > nres2 {
                continue;
            }
            let mut x = 0.0f32;
            for &i0 in &aligned {
                let j0 = preali1[i0] as usize;
                if i != i0 && j != j0 {
                    x += scorefun_w(d1[[i0, i]], d2[[j0, j]], wght);
                }
            }
            score_map.insert((i, j), x);
        }
    }

    // Accept bands +-3 around aligned segments
    for i in 1..=nres1 {
        let k = preali1[i] as i32;
        if k > 0 {
            let j_lo = 1.max(k - 3) as usize;
            let j_hi = ((nres1 as i32 - 3).min(k + 3)) as usize;
            for j in j_lo..=j_hi {
                score_map.insert((i, j), 1.0f32);
            }
        }
    }

    // Count tetrapeptides (4 consecutive positive scores)
    let mut tetrapool = Vec::new();
    if nres1 < 4 {
        return tetrapool;
    }
    for i in 1..=(nres1 - 3) {
        for j_i32 in (left[i] + 1)..rite[i] {
            let j = j_i32 as usize;
            if j < 1 || j > nres2 - 3 {
                continue;
            }
            let mut x = 0.0f32;
            for l in 0..4 {
                x += score_map.get(&(i + l, j + l)).copied().unwrap_or(0.0);
            }
            if x > 0.0 && tetrapool.len() < MAXPAIR {
                tetrapool.push((i, j));
            }
        }
    }

    tetrapool
}

// ---------------------------------------------------------------------------
// Lean Monte Carlo state
// ---------------------------------------------------------------------------

struct LeanState {
    ali1: Vec<i16>,         // nres1+2, 1-based alignment
    ali2: Vec<i16>,         // nres2+2, 1-based reverse
    fragali1: Vec<i16>,     // gene->j mapping
    fragali2: Vec<i16>,     // gene reverse
    cnt: Vec<i32>,          // gene overlap count
    prev1: Vec<i32>,        // linked list prev
    next1: Vec<i32>,        // linked list next
    nextres: Vec<i32>,      // aligned residue chain
    next2: Vec<i32>,        // cd2 linked list
    ltest: Vec<bool>,       // (nres1+2)*(nres2+2) flat
    rescore: Vec<f32>,      // (nres1+2)*(nres2+2) flat
    dscore: Vec<f32>,       // (nres1+2)*(nres2+2) flat
    ncand: usize,
    candij_i: Vec<i32>,     // candidate i indices
    candij_j: Vec<i32>,     // candidate j indices
    nres1: usize,
    nres2: usize,
    stride: usize,          // nres2+2 for flat indexing
}

impl LeanState {
    fn new(nres1: usize, nres2: usize) -> Self {
        let sz1 = nres1 + 2;
        let sz2 = nres2 + 2;
        let flat_sz = sz1 * sz2;
        let max_cand = MAXPAIR * 4 + 1;

        LeanState {
            ali1: vec![0i16; sz1],
            ali2: vec![0i16; sz2],
            fragali1: vec![0i16; sz1],
            fragali2: vec![0i16; sz2],
            cnt: vec![0i32; sz1],
            prev1: vec![0i32; sz1],
            next1: vec![nres1 as i32 + 1; sz1],
            nextres: vec![nres1 as i32 + 1; sz1],
            next2: vec![nres2 as i32 + 1; sz2],
            ltest: vec![false; flat_sz],
            rescore: vec![0.0f32; flat_sz],
            dscore: vec![0.0f32; flat_sz],
            ncand: 0,
            candij_i: vec![0i32; max_cand],
            candij_j: vec![0i32; max_cand],
            nres1,
            nres2,
            stride: sz2,
        }
    }

    #[inline]
    fn flat_idx(&self, i: usize, j: usize) -> usize {
        i * self.stride + j
    }
}

// ---------------------------------------------------------------------------
// Lean MC subroutines
// ---------------------------------------------------------------------------

/// Marginal score of adding pair (i,j).
fn addscore(
    i: usize,
    j: usize,
    state: &LeanState,
    d1: &Array2<i16>,
    d2: &Array2<i16>,
    wght: &[f32; 101],
) -> f32 {
    let zero = 0i16;
    let mut dx = scorefun_w(zero, zero, wght);
    let iold = state.ali2[j] as usize;
    let jold = state.ali1[i] as i32;
    let mut k = 0usize;
    while (state.nextres[k] as usize) <= state.nres1 {
        k = state.nextres[k] as usize;
        let x_val = state.ali1[k] as i32;
        if k != iold && x_val != jold {
            let s = scorefun_w(d1[[i, k]], d2[[j, x_val as usize]], wght);
            // Fortran: dx=dx+s+s → (dx+s)+s
            dx = (dx + s) + s;
        }
    }
    dx
}

/// Calculate rescore on the fly.
fn getrescore(
    i1: usize,
    i2: usize,
    state: &LeanState,
    d1: &Array2<i16>,
    d2: &Array2<i16>,
    wght: &[f32; 101],
) -> f32 {
    let zero = 0i16;
    let mut dx = scorefun_w(zero, zero, wght);
    let mut k = 0usize;
    while (state.nextres[k] as usize) <= state.nres1 {
        k = state.nextres[k] as usize;
        let x_val = state.ali1[k] as usize;
        if k != i1 && x_val != i2 {
            let s = scorefun_w(d1[[i1, k]], d2[[i2, x_val]], wght);
            dx = (dx + s) + s;
        }
    }
    dx
}

/// Add candidates to candidate list.
fn appendcand_lean(
    i: usize,
    j: usize,
    blocksize: usize,
    state: &mut LeanState,
    d1: &Array2<i16>,
    d2: &Array2<i16>,
    wght: &[f32; 101],
    lsave: bool,
) {
    if !lsave {
        return;
    }
    for l in 0..blocksize {
        let j1 = (j + l) as usize; // j is always positive in DALICON
        let idx = state.flat_idx(i + l, j1);
        if !state.ltest[idx] {
            state.ncand += 1;
            let nc = state.ncand;
            state.candij_i[nc] = (i + l) as i32;
            state.candij_j[nc] = j1 as i32;
            state.ltest[idx] = true;

            if state.ali1[i + l] != j1 as i16 {
                let ds = addscore(i + l, j1, state, d1, d2, wght);
                let idx2 = state.flat_idx(i + l, j1);
                state.dscore[idx2] = ds;
            }
        }
    }
}

/// Topology check: find genes violating sequential ordering.
fn checktopo(
    i: usize,
    j: usize,
    blocksize: usize,
    state: &LeanState,
) -> Vec<usize> {
    let mut genedel = Vec::new();
    let nres1 = state.nres1;

    // Forward check
    let mut k = (nres1 + 1).min(i + blocksize - 1);
    while (state.next1[k] as usize) <= nres1 {
        k = state.next1[k] as usize;
        let l = state.fragali1[k] as i32;
        if l <= j as i32 - blocksize as i32 {
            genedel.push(k);
        } else {
            break;
        }
    }

    // Backward check
    let mut k = if i + 1 >= blocksize { i + 1 - blocksize } else { 0 };
    while state.prev1[k] > 0 {
        k = state.prev1[k] as usize;
        let l = state.fragali1[k] as i32;
        if l >= j as i32 + blocksize as i32 {
            genedel.push(k);
        } else {
            break;
        }
    }

    genedel
}

/// Find overlapping genes.
fn geneoverlap(
    i: usize,
    j: usize,
    blocksize: usize,
    state: &LeanState,
    mut genedel: Vec<usize>,
) -> Vec<usize> {
    let nres1 = state.nres1;
    let nres2 = state.nres2;

    // i direction
    let start = if i >= blocksize { i - blocksize } else { 0 };
    let mut k = state.next1[start] as usize;
    while k <= nres1.min(i + blocksize - 1) {
        let l = state.fragali1[k] as i32;
        if (k as i32 - l) != (i as i32 - j as i32) {
            genedel.push(k);
        }
        k = state.next1[k] as usize;
    }

    // j direction (parallel only in DALICON)
    if j > 0 {
        // parallel-parallel
        let start_j = if j >= blocksize { j - blocksize } else { 0 };
        let mut l = state.next2[start_j] as usize;
        while l <= nres2.min(j + blocksize - 1) {
            let k2 = state.fragali2[l] as i32;
            if k2 > 0
                && (k2 - l as i32) != (i as i32 - j as i32)
                && (k2 >= i as i32 + blocksize as i32 || k2 <= i as i32 - blocksize as i32)
            {
                genedel.push(k2 as usize);
            }
            l = state.next2[l] as usize;
        }

        // parallel-antiparallel: always clash
        let start_j2 = if j >= 1 { j - 1 } else { 0 };
        let mut l = state.next2[start_j2] as usize;
        while l <= nres2.min(j + 2 * blocksize - 1) {
            let k2 = state.fragali2[l] as i32;
            if k2 < 0
                && (-k2 >= i as i32 + blocksize as i32
                    || -k2 <= i as i32 - blocksize as i32)
            {
                genedel.push((-k2) as usize);
            }
            l = state.next2[l] as usize;
        }
    }

    genedel
}

/// Find residues that drop to cnt==0 after gene deletion.
fn gowithgenes(
    i1: usize,
    j1: i32,
    genedel: &[usize],
    state: &LeanState,
    blocksize: usize,
    laddition: bool,
) -> Vec<(usize, usize)> {
    // Copy cnt for affected positions
    let mut tmp = std::collections::HashMap::new();
    for &gene_i in genedel {
        for l in 0..blocksize {
            let k = gene_i + l;
            tmp.entry(k).or_insert(state.cnt[k]);
        }
    }

    let mut rem = Vec::new();
    for &gene_i in genedel {
        for l in 0..blocksize {
            let k = gene_i + l;
            *tmp.get_mut(&k).unwrap() -= 1;
            if tmp[&k] == 0 {
                if laddition {
                    let m = k as i32 - i1 as i32;
                    if m >= 0
                        && (m as usize) < blocksize
                        && state.ali1[k] == (j1.abs() + m) as i16
                    {
                        continue;
                    }
                }
                rem.push((k, state.ali1[k] as usize));
            }
        }
    }

    rem
}

/// Test adding gene (i,j): compute score change and affected residues.
fn testaddition(
    i: usize,
    j: usize,
    blocksize: usize,
    state: &LeanState,
    d1: &Array2<i16>,
    d2: &Array2<i16>,
    wght: &[f32; 101],
    ltop: bool,
    lsave: bool,
) -> (f32, Vec<(usize, usize)>, Vec<(usize, usize)>, Vec<usize>) {
    let mut genedel = Vec::new();

    if ltop {
        genedel = checktopo(i, j, blocksize, state);
    }
    genedel = geneoverlap(i, j, blocksize, state, genedel);

    let rem = gowithgenes(i, j as i32, &genedel, state, blocksize, true);

    // New residues
    let mut new_pairs = Vec::new();
    for l in 0..blocksize {
        if state.ali1[i + l] != (j + l) as i16 {
            new_pairs.push((i + l, j + l));
        }
    }

    // Calculate dx (f32)
    let mut dx = 0.0f32;

    for (inew_idx, &(i1, i2)) in new_pairs.iter().enumerate() {
        if lsave {
            dx += state.dscore[state.flat_idx(i1, i2)];
        } else {
            dx += getrescore(i1, i2, state, d1, d2, wght);
        }

        // new × new
        for k in 0..inew_idx {
            let s = scorefun_w(d1[[i1, new_pairs[k].0]], d2[[i2, new_pairs[k].1]], wght);
            dx = (dx + s) + s;
        }

        // new × rem (subtract)
        for &(remi, remj) in &rem {
            if i1 != remi && i2 != remj {
                let s = scorefun_w(d1[[i1, remi]], d2[[i2, remj]], wght);
                dx = (dx - s) - s;
            }
        }
    }

    // rem × rem (subtract rescores, add cross-terms back)
    for (k_idx, &(remi, remj)) in rem.iter().enumerate() {
        if lsave {
            dx -= state.rescore[state.flat_idx(remi, remj)];
        } else {
            dx -= getrescore(remi, remj, state, d1, d2, wght);
        }
        for l_idx in 0..k_idx {
            let s = scorefun_w(d1[[remi, rem[l_idx].0]], d2[[remj, rem[l_idx].1]], wght);
            dx = (dx + s) + s;
        }
    }

    (dx, rem, new_pairs, genedel)
}

/// Test deleting genes: compute score change.
fn testdeletion(
    genedel: &[usize],
    blocksize: usize,
    state: &LeanState,
    d1: &Array2<i16>,
    d2: &Array2<i16>,
    wght: &[f32; 101],
    lsave: bool,
) -> (f32, Vec<(usize, usize)>) {
    let mut all_rem = Vec::new();
    let mut dx = 0.0f32;

    for &gene_i in genedel {
        let j = state.ali1[gene_i] as i32;
        let rem = gowithgenes(gene_i, j, &[gene_i], state, blocksize, false);

        for (l_idx, &(i1, i2)) in rem.iter().enumerate() {
            if lsave {
                dx -= state.rescore[state.flat_idx(i1, i2)];
            } else {
                dx -= getrescore(i1, i2, state, d1, d2, wght);
            }
            for k_idx in 0..l_idx {
                let (j1, j2) = rem[k_idx];
                let s = scorefun_w(d1[[i1, j1]], d2[[i2, j2]], wght);
                dx = (dx + s) + s;
            }
        }

        all_rem.extend(rem);
    }

    (dx, all_rem)
}

/// Execute gene and residue deletions.
fn dodeletions_lean(
    genedel: &[usize],
    rem: &[(usize, usize)],
    state: &mut LeanState,
    blocksize: usize,
    d1: &Array2<i16>,
    d2: &Array2<i16>,
    wght: &[f32; 101],
    lsave: bool,
) -> usize {
    // Delete genes
    for &gene_i in genedel {
        let j = state.fragali1[gene_i] as i32;
        state.fragali1[gene_i] = 0;
        state.fragali2[j.unsigned_abs() as usize] = 0;
        for l in 0..blocksize {
            state.cnt[gene_i + l] -= 1;
        }

        // Update prev1, next1
        let ip = state.prev1[gene_i];
        let in_ = state.next1[gene_i];
        let in_usize = in_ as usize;
        for ix in gene_i..=in_usize.min(state.prev1.len() - 1) {
            state.prev1[ix] = ip;
        }
        for ix in (ip as usize)..=gene_i {
            if ix < state.next1.len() {
                state.next1[ix] = in_;
            }
        }

        // Update next2
        let j1 = j.unsigned_abs() as usize;
        let in_ = state.next2[j1];
        let mut ix = j1 as i32;
        if ix > 0 {
            while ix > 0 && state.next2[(ix - 1) as usize] == j1 as i32 {
                ix -= 1;
                state.next2[ix as usize] = in_;
                if ix <= 0 {
                    break;
                }
            }
        }
    }

    // Delete residues
    let mut ndel = 0;
    let ncand = state.ncand;
    for &(ri, rj) in rem {
        // Update nextres
        let in_ = state.nextres[ri];
        let mut ix = ri as i32;
        while ix > 0 && state.nextres[(ix - 1) as usize] == ri as i32 {
            ix -= 1;
            state.nextres[ix as usize] = in_;
            if ix <= 0 {
                break;
            }
        }

        if lsave {
            // Update dscore/rescore
            for k in 1..=ncand {
                let cix = state.candij_i[k] as usize;
                let cjx = state.candij_j[k] as usize;
                let mut ddx = 0.0f32;
                if ri != cix && rj != cjx {
                    let s = -scorefun_w(d1[[cix, ri]], d2[[cjx, rj]], wght);
                    ddx = s + s;
                }
                let idx = state.flat_idx(cix, cjx);
                if state.ali1[cix] == cjx as i16 {
                    state.rescore[idx] += ddx;
                } else {
                    state.dscore[idx] += ddx;
                }
            }
            let idx = state.flat_idx(ri, rj);
            state.dscore[idx] = state.rescore[idx];
            state.rescore[idx] = 0.0f32;
        }

        state.ali1[ri] = 0;
        state.ali2[rj] = 0;
        ndel += 1;
    }

    ndel
}

/// Execute gene and residue additions.
fn doadditions_lean(
    i: usize,
    j: usize,
    new_pairs: &[(usize, usize)],
    state: &mut LeanState,
    blocksize: usize,
    d1: &Array2<i16>,
    d2: &Array2<i16>,
    wght: &[f32; 101],
    lsave: bool,
) -> usize {
    // Add gene
    state.fragali1[i] = j as i16;
    if j > 0 {
        state.fragali2[j] = i as i16;
    }
    for l in 0..blocksize {
        state.cnt[i + l] += 1;
    }

    // Update prev1, next1
    let ip = state.prev1[i] as usize;
    let in_ = state.next1[i] as usize;
    for ix in ip..i {
        state.next1[ix] = i as i32;
    }
    for ix in (i + 1)..=in_.min(state.prev1.len() - 1) {
        state.prev1[ix] = i as i32;
    }

    // Update next2
    let j1 = j;
    let in_ = state.next2[j1];
    let mut ix = j1 as i32;
    if ix > 0 {
        while ix > 0 && state.next2[(ix - 1) as usize] == in_ {
            ix -= 1;
            state.next2[ix as usize] = j1 as i32;
            if ix <= 0 {
                break;
            }
        }
    }

    // Add residues
    let ncand = state.ncand;
    let mut nacc = 0;
    for &(i1, j1_val) in new_pairs {
        // Update nextres
        let in_ = state.nextres[i1];
        let mut ix = i1 as i32;
        while ix > 0 && state.nextres[(ix - 1) as usize] == in_ {
            ix -= 1;
            state.nextres[ix as usize] = i1 as i32;
            if ix <= 0 {
                break;
            }
        }

        if lsave {
            for k in 1..=ncand {
                let cix = state.candij_i[k] as usize;
                let cjx = state.candij_j[k] as usize;
                let mut ddx = 0.0f32;
                if i1 != cix && j1_val != cjx {
                    let s = scorefun_w(d1[[cix, i1]], d2[[cjx, j1_val]], wght);
                    ddx = s + s;
                }
                let idx = state.flat_idx(cix, cjx);
                if state.ali1[cix] == cjx as i16 {
                    state.rescore[idx] += ddx;
                } else {
                    state.dscore[idx] += ddx;
                }
            }
            let idx = state.flat_idx(i1, j1_val);
            state.rescore[idx] = state.dscore[idx];
            state.dscore[idx] = 0.0f32;
        }

        // Handle conflicts
        let iold = state.ali2[j1_val] as usize;
        let jold = state.ali1[i1] as usize;
        if iold != 0 {
            state.ali1[iold] = 0;
        }
        if jold != 0 {
            state.ali2[jold] = 0;
        }
        state.ali1[i1] = j1_val as i16;
        state.ali2[j1_val] = i1 as i16;
        nacc += 1;
    }

    nacc
}

/// Apply changes: deletions then additions.
fn changeconfig_lean(
    genedel: &[usize],
    rem: &[(usize, usize)],
    i: usize,
    j: usize,
    new_pairs: &[(usize, usize)],
    state: &mut LeanState,
    blocksize: usize,
    d1: &Array2<i16>,
    d2: &Array2<i16>,
    wght: &[f32; 101],
    dx: f32,
    totscore: f32,
    bestscore: &mut f32,
    bestali1: &mut [i16],
    outfragali: &mut [i16],
    lsave: bool,
) -> (usize, usize, f32) {
    let ndel = dodeletions_lean(genedel, rem, state, blocksize, d1, d2, wght, lsave);
    let mut nacc = 0;
    if i != 0 {
        nacc = doadditions_lean(i, j, new_pairs, state, blocksize, d1, d2, wght, lsave);
    }

    let new_totscore = totscore + dx;

    if new_totscore > *bestscore {
        // Save best fragali and ali1
        for k in 1..=state.nres1 {
            outfragali[k] = state.fragali1[k];
            bestali1[k] = state.ali1[k];
        }
        *bestscore = new_totscore;
    }

    (ndel, nacc, new_totscore)
}

/// Check if gene at position i has cnt==1 at either end.
#[inline]
fn cnt_at_end(cnt: &[i32], i: usize, blocksize: usize) -> bool {
    cnt[i] == 1 || cnt[i + blocksize - 1] == 1
}

/// Non-recursive quicksort of index array d by values in ve.
/// Sorts d[lo..=hi] (1-based).
fn index_qsort(d: &mut [i32], ve: &[i32], lo: usize, hi: usize) {
    let mut stack = Vec::new();
    stack.push((lo, hi));

    while let Some((mut top, mut bottom)) = stack.pop() {
        while top < bottom {
            let mut upper = top;
            let mut lower = bottom;
            let save = d[upper];
            while upper != lower {
                while upper < lower && ve[save as usize] <= ve[d[lower] as usize] {
                    lower -= 1;
                }
                if upper != lower {
                    d[upper] = d[lower];
                }
                while upper < lower && ve[save as usize] >= ve[d[upper] as usize] {
                    upper += 1;
                }
                if upper != lower {
                    d[lower] = d[upper];
                }
            }
            d[upper] = save;
            let p = upper;

            if (p - top) > (bottom - p) {
                stack.push((top, if p > 0 { p - 1 } else { 0 }));
                top = p + 1;
            } else {
                stack.push((p + 1, bottom));
                bottom = if p > 0 { p - 1 } else { 0 };
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Main lean_mc
// ---------------------------------------------------------------------------

/// Main lean Monte Carlo optimization.
///
/// Returns (score, bestali1) where score = bestscore / 100.
fn lean_mc(
    imax0: usize,
    blocksize: usize,
    nres1: usize,
    nres2: usize,
    d1: &Array2<i16>,
    d2: &Array2<i16>,
    lsave: bool,
    preali1: &[i16],
    _itrimx: usize,
    ntetra: usize,
    tetrapool: &[(usize, usize)],
    wght: &[f32; 101],
) -> (f32, Vec<i16>) {
    let imax = imax0;
    let sz1 = nres1 + 2;
    let mut bestali1 = vec![0i16; sz1];
    let mut outfragali = vec![0i16; sz1];
    let mut bestscore = 0.0f32;

    let mut state = LeanState::new(nres1, nres2);

    // Load candidates
    for ix in 0..ntetra {
        let (i, j) = tetrapool[ix];
        appendcand_lean(i, j, blocksize, &mut state, d1, d2, wght, lsave);
    }

    // Convert preali1 to initfragali
    let mut initfragali = vec![0i16; sz1];
    let preali_len = preali1.len();
    for i in 1..nres1.min(preali_len).saturating_sub(2) {
        let j = preali1[i] as i32;
        if j > 0 && j <= nres2 as i32 - 3 {
            if i > 1 && preali1[i - 1] == 0 {
                initfragali[i] = j as i16;
            }
            if i + 3 < preali_len && preali1[i + 3] == (j + 3) as i16 {
                initfragali[i] = j as i16;
            }
        }
    }

    // Load initial alignment
    let mut totscore = 0.0f32;
    for i in 1..=nres1 {
        let j = initfragali[i] as i32;
        if j != 0 {
            let (dx, rem, new_pairs, genedel) =
                testaddition(i, j as usize, blocksize, &state, d1, d2, wght, true, lsave);
            bestscore = 0.0f32; // Force accept (aconitase fix)
            let (_, _, ts) = changeconfig_lean(
                &genedel, &rem, i, j as usize, &new_pairs,
                &mut state, blocksize, d1, d2, wght,
                dx, totscore, &mut bestscore, &mut bestali1, &mut outfragali, lsave,
            );
            totscore = ts;
        }
    }

    // Monte Carlo steps
    let mut imp = 0usize;
    let mut istep = 0usize;
    let mut impscore = 0.0f32;
    let mut alfa = 50.0f32;
    let mut delta = 0.0f32;

    let mut di = vec![0i32; ntetra + 1];
    let mut ve = vec![0i32; ntetra + 1];

    while istep < imax {
        istep += 1;

        // Heat shock
        if istep - imp == 1 {
            alfa = 50.0f32;
            delta = 0.0f32;
        } else if istep - imp == 2 {
            alfa = 1.0f32 / 10000.0f32.max(totscore).sqrt();
            delta = alfa / 1000.0f32;
        }

        let mut ndel_total = 0usize;
        let mut nacc_total = 0usize;
        alfa += delta;

        // Adding phase (WARM) — randomized order
        for ix in 1..=ntetra {
            di[ix] = ix as i32;
            ve[ix] = nint((100.0f32 * RNG_CONSTANT) as f64);
        }
        if ntetra > 0 {
            index_qsort(&mut di, &ve, 1, ntetra);
        }

        for q in 1..=ntetra {
            let ix = di[q] as usize;
            let (i, j) = tetrapool[ix - 1];
            if state.fragali1[i] != j as i16 {
                let (dx, rem, new_pairs, genedel) =
                    testaddition(i, j, blocksize, &state, d1, d2, wght, true, lsave);
                let p = getp(alfa * dx);
                if RNG_CONSTANT < p {
                    let (ndel, nacc, ts) = changeconfig_lean(
                        &genedel, &rem, i, j, &new_pairs,
                        &mut state, blocksize, d1, d2, wght,
                        dx, totscore, &mut bestscore, &mut bestali1, &mut outfragali, lsave,
                    );
                    totscore = ts;
                    ndel_total += ndel;
                    nacc_total += nacc;
                    if totscore > impscore {
                        imp = istep;
                        impscore = totscore;
                    }
                }
            }
        }

        // Trimming phase (COLD)
        {
            let mut i = 0usize;
            while (state.next1[i] as usize) <= nres1 {
                i = state.next1[i] as usize;
                if cnt_at_end(&state.cnt, i, blocksize) {
                    let genedel_trim = vec![i];
                    let (dx_trim, rem_trim) =
                        testdeletion(&genedel_trim, blocksize, &state, d1, d2, wght, lsave);
                    if dx_trim > 0.0f32 {
                        let (ndel, nacc, ts) = changeconfig_lean(
                            &genedel_trim, &rem_trim, 0, 0, &[],
                            &mut state, blocksize, d1, d2, wght,
                            dx_trim, totscore, &mut bestscore, &mut bestali1, &mut outfragali, lsave,
                        );
                        totscore = ts;
                        ndel_total += ndel;
                        nacc_total += nacc;
                        if totscore > impscore + 1.0f32 {
                            imp = istep;
                            impscore = totscore;
                        }
                    }
                }
            }
        }

        // Early return
        let ncand = state.ncand;
        if nacc_total > ncand || ndel_total > ncand {
            return (bestscore / 100.0, bestali1);
        }

        // Quit if no improvement for 20 steps
        if istep >= imp + 20 {
            break;
        }

        // Accelerate cooling
        if ntetra > 0 && (nacc_total as f64 / ntetra as f64) > 0.1 {
            alfa = 10.0f32 * alfa;
        }
    }

    (bestscore / 100.0, bestali1)
}

// ---------------------------------------------------------------------------
// Output: convert alignment to blocks
// ---------------------------------------------------------------------------

/// Extract block assignments from alignment.
///
/// getblocks with blocksize=1. Uses two separate `if` (not `else if`).
fn getblocks(ali1: &[i16], nres1: usize) -> (Vec<i32>, usize) {
    let mut block = vec![0i32; nres1 + 2];
    let mut nblock = 0usize;

    if nres1 >= 1 && ali1[1] != 0 {
        nblock = 1;
        block[1] = 1;
    }

    for i in 2..=nres1 {
        if ali1[i] != 0 {
            // Two SEPARATE if checks — both can increment nblock
            if ali1[i] == 1 {
                nblock += 1;
            }
            if ali1[i] != ali1[i - 1] + 1 && ali1[i] != ali1[i - 1] - 1 {
                nblock += 1;
            }
            block[i] = nblock as i32;
        }
    }

    (block, nblock)
}

/// Filter: set score to 0 if too few blocks AND too few residues.
fn filter_score(alilen: usize, nblock: usize, nres1: usize, nres2: usize) -> bool {
    let minlen = (nres1 / 2).min(nres2 / 2).min(30);
    !(alilen < minlen && nblock < 4)
}

/// Convert alignment to output blocks.
///
/// Returns (nblock, blocks) or None if filtered out.
fn produce_output(
    nres1: usize,
    nres2: usize,
    ali1: &[i16],
    score: f32,
) -> Option<(usize, Vec<AlignmentBlock>)> {
    // getali0 with blocksize=1 is identity
    let a1 = ali1;

    let (block, nblock) = getblocks(ali1, nres1);

    // Compute alignment length
    let lenali = (1..=nres1).filter(|&j| a1[j] != 0).count();

    if score <= 0.0 || !filter_score(lenali, nblock, nres1, nres2) {
        return None;
    }

    // Build block boundaries
    let mut ib1 = vec![0usize; nblock + 1];
    let mut ib2 = vec![0usize; nblock + 1];
    for j in 1..=nres1 {
        let b = block[j] as usize;
        if b > 0 {
            if ib1[b] == 0 {
                ib1[b] = j;
            }
            ib2[b] = j;
        }
    }

    let mut blocks = Vec::new();
    for b in 1..=nblock {
        blocks.push(AlignmentBlock {
            l1: ib1[b] as u32,
            r1: ib2[b] as u32,
            l2: a1[ib1[b]] as u32,
            r2: a1[ib2[b]] as u32,
        });
    }

    Some((nblock, blocks))
}

// ---------------------------------------------------------------------------
// Prealignment conversion
// ---------------------------------------------------------------------------

/// Convert block ranges to per-residue prealignment array.
///
/// Translated from comparemodules.f dowork_dalicon prealignment setup.
fn convert_prealignment(
    nblock: usize,
    values: &[i32],
    nres1: usize,
) -> (Vec<i16>, usize) {
    let mut preali1 = vec![0i16; nres1 + 1];
    let mut nx = 0usize;

    // Fortran: do i=1,npr,2 where npr=nblock
    //   k=pr(npr*2+i); do j=pr(i),pr(i+1); preali1(j)=k; k=k+1; end do
    // values[0..2*nblock] = cd1 ranges, values[2*nblock..4*nblock] = cd2 ranges
    let npr = nblock;
    for i_fort in (0..npr).step_by(2) {
        let start = values[i_fort] as usize;
        let end = values[i_fort + 1] as usize;
        let mut k = values[npr * 2 + i_fort]; // cd2 start for this range
        for j in start..=end {
            if j >= 1 && j <= nres1 {
                preali1[j] = k as i16;
                k += 1;
                nx += 1;
            }
        }
    }

    (preali1, nx)
}

// ---------------------------------------------------------------------------
// Public interface
// ---------------------------------------------------------------------------

/// Mutable state for cd1 that persists across comparisons.
///
/// In Fortran, ca1 is a module-level variable modified by transrotate in-place.
/// We make this explicit.
pub struct DaliconCd1State {
    pub code: String,
    pub nres: usize,
    pub ca: Array2<f64>,
    pub dist: DaliconDistMat,
}

/// DALICON input record parsed from text file.
pub struct DaliconRecord {
    pub cd1: String,
    pub cd2: String,
    pub nblock: usize,
    pub values: Vec<i32>,
}

/// Parse DALICON input file (grouped by cd1 with END delimiters).
pub fn parse_dalicon_input(content: &str) -> Vec<DaliconRecord> {
    let lines: Vec<&str> = content.lines().collect();
    let mut records = Vec::new();
    let mut i = 0;
    let mut cd1: Option<String> = None;

    while i < lines.len() {
        let line = lines[i].trim();

        if line == "END" {
            i += 1;
            // Skip blank lines
            while i < lines.len() && lines[i].trim().is_empty() {
                i += 1;
            }
            if i < lines.len() && lines[i].trim() != "END" {
                cd1 = Some(lines[i].trim().to_string());
                i += 1;
            }
            continue;
        }

        if line.is_empty() {
            i += 1;
            continue;
        }

        if cd1.is_none() {
            cd1 = Some(line.to_string());
            i += 1;
            continue;
        }

        // Try parsing as cd2 name (not all digits)
        let is_numeric = line
            .replace('-', "")
            .replace('.', "")
            .chars()
            .all(|c| c.is_ascii_digit() || c.is_whitespace());

        if !is_numeric {
            let cd2 = line.trim_end_matches('*').to_string();
            i += 1;
            // Skip blank lines
            while i < lines.len() && lines[i].trim().is_empty() {
                i += 1;
            }
            if i >= lines.len() {
                break;
            }
            let nblock: usize = lines[i].trim().parse().unwrap_or(0);
            i += 1;

            // Read 4*nblock integers
            let mut values = Vec::new();
            while values.len() < nblock * 4 && i < lines.len() {
                let vline = lines[i].trim();
                if vline.is_empty() || vline == "END" {
                    break;
                }
                for tok in vline.split_whitespace() {
                    if let Ok(v) = tok.parse::<i32>() {
                        values.push(v);
                    }
                }
                i += 1;
            }

            records.push(DaliconRecord {
                cd1: cd1.clone().unwrap(),
                cd2,
                nblock,
                values,
            });
        } else {
            i += 1;
        }
    }

    records
}

/// Run DALICON on a single pair.
///
/// Returns output blocks or None.
pub fn dowork_dalicon(
    cd1: &str,
    cd2: &str,
    store: &ProteinStore,
    nblock: usize,
    values: &[i32],
    cd1_state: &mut Option<DaliconCd1State>,
    lfitz: bool,
) -> Option<(usize, Vec<AlignmentBlock>)> {
    let wght = gagaweights();

    // Setup cd1
    let need_new_cd1 = match cd1_state {
        Some(ref s) => s.code != cd1,
        None => true,
    };

    if need_new_cd1 {
        let prot1 = store.get_protein(cd1).ok()?;
        let ca1 = prot1.ca.clone();
        let d1 = getgagadist(&ca1, prot1.nres);
        *cd1_state = Some(DaliconCd1State {
            code: cd1.to_string(),
            nres: prot1.nres,
            ca: ca1,
            dist: d1,
        });
    }

    let cd1_s = cd1_state.as_mut().unwrap();
    let nres1 = cd1_s.nres;

    // Setup cd2
    let prot2 = store.get_protein(cd2).ok()?;
    let nres2 = prot2.nres;
    let ca2 = &prot2.ca;
    let d2 = getgagadist(ca2, nres2);

    // Convert prealignment
    let (mut preali1, nx) = convert_prealignment(nblock, values, nres1);

    // Run fitz to extend short prealignment
    if lfitz {
        let mut ali_0based = vec![-1i32; nres1];
        for i in 1..=nres1 {
            if preali1[i] > 0 {
                ali_0based[i - 1] = preali1[i] as i32 - 1;
            }
        }

        // Initial superposition from prealignment (Python: getut + transrotate)
        {
            let pairs: Vec<(usize, usize)> = (0..nres1)
                .filter(|&i| ali_0based[i] >= 0)
                .map(|i| (i, ali_0based[i] as usize))
                .collect();
            let lali = pairs.len();
            if lali > 0 {
                let mut ux = Array2::zeros((3, lali));
                let mut uy = Array2::zeros((3, lali));
                let w = vec![1.0; lali];
                for (k, &(i_idx, j_idx)) in pairs.iter().enumerate() {
                    for dim in 0..3 {
                        ux[[dim, k]] = cd1_s.ca[[dim, i_idx]];
                        uy[[dim, k]] = ca2[[dim, j_idx]];
                    }
                }
                let result = u3b(&w, &ux, &uy, lali, 1);
                transrotate(&mut cd1_s.ca, &result.u, &result.t);
            }
        }

        // fitz modifies ca1 in-place (via transrotate)
        let _fitz_result = fitz(&mut cd1_s.ca, ca2, FITZRCUT, FITZMAXITER);

        // Overwrite preali1 if fitz result is longer
        let nali_fitz = _fitz_result.lali;
        if nali_fitz > nx {
            for i in 0..nres1 {
                if _fitz_result.ali[i] >= 0 {
                    preali1[i + 1] = (_fitz_result.ali[i] + 1) as i16;
                } else {
                    preali1[i + 1] = 0;
                }
            }
        }
    }

    // Testi: find tetrapeptide candidates
    let tetrapool = testi(nres1, nres2, &preali1, &cd1_s.dist.data, &d2.data, wght);
    let ntetra = tetrapool.len();

    // Initialize
    let mut score = 0.0f32;
    let mut ali1 = vec![0i16; MAXRES + 1];

    if ntetra > 0 {
        let (s, a) = lean_mc(
            IMAX0, BLOCKSIZE, nres1, nres2,
            &cd1_s.dist.data, &d2.data, true,
            &preali1, ITRIM2, ntetra, &tetrapool, wght,
        );
        score = s;
        let copy_len = ali1.len().min(a.len());
        ali1[..copy_len].copy_from_slice(&a[..copy_len]);
    }

    // Iterate until improvement < refitol
    let mut oldscore = score;
    let mut ds = oldscore;
    let mut iter_count = 0;

    while ds > REFITOL && iter_count < MAXITER {
        iter_count += 1;
        for i in 1..=nres1 {
            preali1[i] = ali1[i];
        }

        let tetrapool = testi(nres1, nres2, &preali1, &cd1_s.dist.data, &d2.data, wght);
        let ntetra = tetrapool.len();

        let mut score_new = 0.0f32;
        let mut ali1_new = vec![0i16; MAXRES + 1];

        if ntetra > 0 {
            let (s, a) = lean_mc(
                IMAX0, BLOCKSIZE, nres1, nres2,
                &cd1_s.dist.data, &d2.data, true,
                &preali1, ITRIM2, ntetra, &tetrapool, wght,
            );
            score_new = s;
            let copy_len = ali1_new.len().min(a.len());
            ali1_new[..copy_len].copy_from_slice(&a[..copy_len]);
        }

        ali1 = ali1_new;
        score = score_new;
        ds = score - oldscore;
        oldscore = score;
    }

    produce_output(nres1, nres2, &ali1, score)
}

/// Run DALICON on all pairs from parsed input records.
///
/// Returns results as Vec<(cd1, cd2, nblock, blocks)>.
pub fn run_dalicon(
    records: &[DaliconRecord],
    store: &ProteinStore,
) -> Vec<(String, String, usize, Vec<AlignmentBlock>)> {
    let mut results = Vec::new();
    let mut cd1_state: Option<DaliconCd1State> = None;

    for rec in records {
        if let Some((nblock, blocks)) = dowork_dalicon(
            &rec.cd1, &rec.cd2, store,
            rec.nblock, &rec.values,
            &mut cd1_state, true,
        ) {
            results.push((rec.cd1.clone(), rec.cd2.clone(), nblock, blocks));
        }
    }

    results
}
