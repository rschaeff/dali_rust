//! SOAP: structural alignment for proteins with fewer than 3 SSEs.
//!
//! Port of DaliLite's `soap4.f` and `dowork_soap` (comparemodules.f).
//! dali.pl runs SOAP in place of WOLF (and without PARSI) when the query has
//! nsse < MINSSE = 3; its WOLFITZ output goes through the same
//! DP → DALICON → DP chain as WOLF's.
//!
//! Arithmetic is REAL (f32) as in Fortran; tables are integer*2. Arrays are
//! kept 1-based (index 0 unused) so the code reads against the Fortran.

use ndarray::Array2;

use crate::numerics::fitz::fitz;
use crate::numerics::kabsch::{transrotate_f32, u3b_f32};
use crate::numerics::nint_f32;
use crate::types::{compress_blocks, AlignmentBlock, Protein};
use crate::wolf::WolfResult;

/// serialcompare.f: rcut, maxiter passed to dowork_soap.
const RCUT_FITZ: f32 = 4.0;
const MAXITER: usize = 20;

/// Query state carried across cd2: Fortran loads ca1 only when cd1 changes
/// (`oldcd1`) and transforms it in place for every domain pair.
pub struct SoapQuery {
    pub code: String,
    nres: usize,
    domains: Vec<Vec<usize>>, // residues (1-based) of each (+,*)-unit
    ca: Array2<f32>,          // (3, nres+1), column 0 unused
}

impl SoapQuery {
    pub fn new(p: &Protein) -> Self {
        SoapQuery {
            code: p.code.clone(),
            nres: p.nres,
            domains: unit_residues(p),
            ca: ca_f32(p),
        }
    }
}

fn ca_f32(p: &Protein) -> Array2<f32> {
    let mut ca = Array2::<f32>::zeros((3, p.nres + 1));
    for i in 0..p.nres {
        for k in 0..3 {
            ca[[k, i + 1]] = p.ca[[k, i]] as f32;
        }
    }
    ca
}

/// setup: keep node 1 and the (+,*)-units; list each unit's residues.
fn unit_residues(p: &Protein) -> Vec<Vec<usize>> {
    p.filtered_domains()
        .iter()
        .map(|d| {
            d.segments
                .iter()
                .flat_map(|&(a, b)| (a as usize)..=(b as usize))
                .collect()
        })
        .collect()
}

/// initequal_fast_domains. `domres2` is the caller's buffer: Fortran indexes
/// it by b without checking b <= domnres2, so a stale entry may be read.
fn initequal_fast_domains(
    nres1: usize,
    nres2: usize,
    ali: &mut [i32],
    domres1: &[usize],
    domnres1: usize,
    domres2: &[usize],
    domnres2: usize,
) {
    for i in 1..=nres1 {
        ali[i] = 0;
    }
    let r = if domnres1 < domnres2 {
        domnres1 as f32 / domnres2 as f32
    } else {
        domnres2 as f32 / domnres1 as f32
    };
    let mut fa = 1.0f32;
    let mut fb = 1.0f32;
    let mut a: i64 = 0;
    while a < domnres1 as i64 {
        if domnres1 < domnres2 {
            fa += r;
            fb += 1.0;
        } else {
            fa += 1.0;
            fb += r;
        }
        a = fa as i64;
        let b = fb as i64;
        if a > 0 && b > 0 && a <= nres1 as i64 && b <= nres2 as i64 {
            let ia = domres1[a as usize];
            if ali[ia] == 0 {
                ali[ia] = domres2[b as usize] as i32;
            }
        }
    }
}

/// getut_w with table0 = 0 (as left by initequal_fast_domains): every pair
/// has weight 1. Returns (u, t).
fn getut_w(nres1: usize, nres2: usize, ali: &[i32], x: &Array2<f32>, y: &Array2<f32>) -> ([[f32; 3]; 3], [f32; 3]) {
    let mut ux: Vec<[f32; 3]> = Vec::new();
    let mut uy: Vec<[f32; 3]> = Vec::new();
    for a in 1..=nres1 {
        let b = (ali[a].max(0) as usize).min(nres2);
        if b == 0 {
            continue;
        }
        ux.push([x[[0, a]], x[[1, a]], x[[2, a]]]);
        uy.push([y[[0, b]], y[[1, b]], y[[2, b]]]);
    }
    let w = vec![1.0f32; ux.len()];
    let r = u3b_f32(&w, &to_mat(&ux), &to_mat(&uy), ux.len(), 1);
    (r.u, r.t)
}

fn to_mat(v: &[[f32; 3]]) -> Array2<f32> {
    let mut m = Array2::<f32>::zeros((3, v.len()));
    for (i, p) in v.iter().enumerate() {
        for k in 0..3 {
            m[[k, i]] = p[k];
        }
    }
    m
}

/// fill1_fast: penalty nint(100*d^2/rcutoff^2), 100 beyond the cutoff.
fn fill1_fast(nres1: usize, nres2: usize, ca1: &Array2<f32>, ca2: &Array2<f32>, table0: &mut Array2<i16>, rcutoff: f32) {
    let r2 = rcutoff * rcutoff;
    for i in 1..=nres1 {
        for j in 1..=nres2 {
            let d1 = ca1[[0, i]] - ca2[[0, j]];
            let d2 = ca1[[1, i]] - ca2[[1, j]];
            let d3 = ca1[[2, i]] - ca2[[2, j]];
            let r = d1 * d1 + d2 * d2 + d3 * d3;
            table0[[i, j]] = if r >= r2 { 100 } else { nint_f32(100.0 * r / r2) as i16 };
        }
    }
}

/// soap_nw_fast: minimum-penalty path from (1,1) to (m,n), filled by
/// anti-diagonals. Returns ali(1..m+n) of moves (1: x, 2: y; a diagonal
/// step is written as 2 then 1), unused leading entries 0.
fn soap_nw_fast(m: usize, n: usize, table0: &Array2<i16>, trace: &mut Array2<i16>) -> Vec<i32> {
    let mut shift = vec![[0i32; 3]; m.max(2) + 1]; // shift[k][diag]
    let (mut diag0, mut diag1, mut diag2) = (2usize, 1usize, 0usize);
    trace[[1, 1]] = 0;
    trace[[2, 1]] = 1;
    trace[[1, 2]] = 2;
    shift[0][diag2] = table0[[1, 1]] as i32;
    shift[0][diag1] = table0[[2, 1]] as i32 + shift[0][diag2];
    shift[1][diag1] = table0[[1, 2]] as i32 + shift[0][diag2];
    for y0 in 3..=(m + n - 1) {
        let k0 = y0.saturating_sub(n); // x == 1
        let (x, y) = (1 + k0, y0 - k0);
        trace[[x, y]] = 2;
        shift[k0][diag0] = shift[k0][diag1] + table0[[x, y]] as i32;
        let k1 = (m - 1).min(y0 - 1); // y == 1
        let (x, y) = (1 + k1, y0 - k1);
        trace[[x, y]] = 1;
        shift[k1][diag0] = shift[k1 - 1][diag1] + table0[[x, y]] as i32;
        for k in (k0 + 1)..k1 {
            let (x, y) = (1 + k, y0 - k);
            let s = table0[[x, y]] as i32;
            let a = shift[k - 1][diag1];
            let b = shift[k][diag1];
            let c = shift[k - 1][diag2];
            if a <= b && a <= c {
                trace[[x, y]] = 1;
                shift[k][diag0] = a + s;
            } else if b <= a && b <= c {
                trace[[x, y]] = 2;
                shift[k][diag0] = b + s;
            } else {
                trace[[x, y]] = 3;
                shift[k][diag0] = c + s;
            }
        }
        diag0 = (diag0 + 1) % 3;
        diag1 = (diag1 + 1) % 3;
        diag2 = (diag2 + 1) % 3;
    }

    let mut ali = vec![0i32; m + n + 1];
    let (mut x, mut y, mut v) = (m, n, m + n);
    while x > 1 || y > 1 {
        let i = trace[[x, y]];
        ali[v] = i as i32;
        v -= 1;
        match i {
            1 => x -= 1,
            2 => y -= 1,
            3 => {
                x -= 1;
                y -= 1;
                ali[v + 1] = 1; // patch 3 -> (1,2)
                ali[v] = 2;
                v -= 1;
            }
            _ => return ali,
        }
    }
    ali
}

/// soap_getut: weighted superposition of the path's pairs. Returns
/// (lali, rms, totw); totw adds each weight three times, as in Fortran.
fn soap_getut(nx: usize, ali: &[i32], x: &Array2<f32>, y: &Array2<f32>, table0: &Array2<i16>) -> (usize, f32, f32) {
    let (mut a, mut b) = (1usize, 1usize);
    let mut totw = 0.0f32;
    let mut ux: Vec<[f32; 3]> = Vec::new();
    let mut uy: Vec<[f32; 3]> = Vec::new();
    let mut w: Vec<f32> = Vec::new();
    for i in 2..=nx.saturating_sub(2) {
        let k = ali[i];
        if k == 0 {
            continue;
        }
        match k {
            1 => a += 1,
            2 => b += 1,
            3 => {
                a += 1;
                b += 1;
            }
            _ => {}
        }
        let s = table0[[a, b]];
        if s == 100 {
            continue;
        }
        let mut px = [0.0f32; 3];
        let mut py = [0.0f32; 3];
        let wl = 1.0f32 - s as f32 / 100.0;
        for j in 0..3 {
            px[j] = if k == 1 { (x[[j, a]] + x[[j, a - 1]]) / 2.0 } else { x[[j, a]] };
            py[j] = if k == 2 { (y[[j, b]] + y[[j, b - 1]]) / 2.0 } else { y[[j, b]] };
            totw += wl;
        }
        ux.push(px);
        uy.push(py);
        w.push(wl);
    }
    let lali = ux.len();
    let ssq = if lali > 2 {
        u3b_f32(&w, &to_mat(&ux), &to_mat(&uy), lali, 1).rms
    } else {
        0.0
    };
    let rms = (ssq / totw.max(1.0)).sqrt();
    (lali, rms, totw)
}

/// dowork_soap for one (cd1, cd2): returns the WOLFITZ alignment.
/// `q.ca` is left transformed, as Fortran's ca1 is.
pub fn soap_compare(q: &mut SoapQuery, cd2: &Protein) -> Option<WolfResult> {
    let nres1 = q.nres;
    let nres2 = cd2.nres;
    if nres1 < 2 || nres2 < 2 {
        return None;
    }
    let ca2 = ca_f32(cd2);
    let domains2 = unit_residues(cd2);
    if q.domains.is_empty() || domains2.is_empty() {
        return None;
    }

    let mut table0 = Array2::<i16>::zeros((nres1 + 1, nres2 + 1));
    let mut trace = Array2::<i16>::zeros((nres1 + 1, nres2 + 1));
    let mut ali = vec![0i32; nres1 + nres2 + 1];
    let mut domres1 = vec![0usize; nres1 + 2];
    let mut domres2 = vec![0usize; nres2 + 2];
    let mut oldw = 0.0f32;
    let mut x1: Option<Array2<f32>> = None;

    for d1 in 0..q.domains.len() {
        let domnres1 = q.domains[d1].len().min(nres1);
        for (i, &r) in q.domains[d1].iter().take(domnres1).enumerate() {
            domres1[i + 1] = r;
        }
        for d2 in &domains2 {
            let domnres2 = d2.len().min(nres2);
            for (i, &r) in d2.iter().take(domnres2).enumerate() {
                domres2[i + 1] = r;
            }
            // init triangulation (table0 := 0)
            initequal_fast_domains(nres1, nres2, &mut ali, &domres1, domnres1, &domres2, domnres2);
            table0.fill(0);
            let (u, t) = getut_w(nres1, nres2, &ali, &q.ca, &ca2);
            transrotate_f32(&mut q.ca, &u, &t);

            // Gaussian penalty; ca1 is not moved inside this loop
            let mut iter = 0;
            let mut r = 30.0f32;
            let mut rms = 99.9f32;
            let mut totw = 0.0f32;
            while iter < 100 && r > 1.0 && rms >= 2.0 {
                iter += 1;
                r *= 0.90;
                fill1_fast(nres1, nres2, &q.ca, &ca2, &mut table0, r);
                let path = soap_nw_fast(nres1, nres2, &table0, &mut trace);
                ali[..path.len()].copy_from_slice(&path);
                let (_lali, new_rms, new_totw) = soap_getut(nres1 + nres2, &ali, &q.ca, &ca2, &table0);
                rms = new_rms;
                totw = new_totw;
            }
            if totw > oldw {
                oldw = totw;
                x1 = Some(q.ca.clone());
            }
        }
    }

    // fitz from the best start; without one Fortran's x1 is undefined
    let x1 = x1?;
    let mut xca = Array2::<f32>::zeros((3, nres1));
    let mut yca = Array2::<f32>::zeros((3, nres2));
    for i in 0..nres1 {
        for k in 0..3 {
            xca[[k, i]] = x1[[k, i + 1]];
        }
    }
    for j in 0..nres2 {
        for k in 0..3 {
            yca[[k, j]] = ca2[[k, j + 1]];
        }
    }
    let fr = fitz(&mut xca, &yca, RCUT_FITZ, MAXITER);
    let blocks: Vec<AlignmentBlock> = (0..nres1)
        .filter(|&i| fr.ali[i] >= 0)
        .map(|i| {
            let j = fr.ali[i] as u32 + 1;
            AlignmentBlock { l1: i as u32 + 1, r1: i as u32 + 1, l2: j, r2: j }
        })
        .collect();
    if blocks.is_empty() {
        return None;
    }
    Some(WolfResult {
        cd1: q.code.clone(),
        cd2: cd2.code.clone(),
        blocks: compress_blocks(&blocks),
    })
}
