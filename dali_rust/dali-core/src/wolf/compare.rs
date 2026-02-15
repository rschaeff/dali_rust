use std::collections::HashMap;
use ndarray::Array2;

use super::geometry::{preparex, twist};
use super::spatial_hash::{SpatialHashGrid, fung};
use crate::types::SseType;

/// Result of SSE pair comparison: best (cd2_i, cd2_j, cd1_u, cd1_v) mapping.
pub struct CompareResult {
    pub iseg: usize,  // cd2 SSE pair first
    pub jseg: usize,  // cd2 SSE pair second
    pub useg: usize,  // cd1 SSE pair first
    pub vseg: usize,  // cd1 SSE pair second
    pub count: usize,  // best vote count
}

/// Find best SSE pair mapping between cd2 (query) and cd1 (in grid).
///
/// For each cd2 SSE pair, transform cd2's other SSEs into canonical frame,
/// look up matching cd1 entries in the spatial hash grid, accumulate votes.
pub fn compare(
    nseg_cd2: usize,
    mid_cd2: &Array2<f64>,
    dir_cd2: &Array2<f64>,
    secstr_cd2: &[SseType],
    nseg_cd1: usize,
    secstr_cd1: &[SseType],
    grid: &SpatialHashGrid,
    rcut: f64,
    neidist_cd2: &Array2<f64>,
) -> Option<CompareResult> {
    let mut protcount: usize = 0;
    let mut best: Option<CompareResult> = None;

    for iseg in 0..nseg_cd2 {
        let atype = secstr_cd2[iseg];
        for jseg in 0..nseg_cd2 {
            if iseg == jseg {
                continue;
            }
            let lest = iseg < jseg; // topology flag for iseg vs jseg ordering

            let r = neidist_cd2[[iseg, jseg]];
            if r >= rcut {
                continue;
            }

            // Vote counts
            let mut count: HashMap<(usize, usize), usize> = HashMap::new();

            // Build canonical frame for cd2 pair (iseg, jseg)
            let mut x = preparex(iseg, jseg, nseg_cd2, mid_cd2, dir_cd2);
            twist(&mut x, 3 + 2 * nseg_cd2);

            // For each other cd2 SSE, look up in cd1's grid
            for kseg in 0..nseg_cd2 {
                if kseg == iseg {
                    continue;
                }
                let less = iseg < kseg; // topology flag for iseg vs kseg ordering

                // Compute midpoint and direction of kseg in canonical frame
                let midx = [
                    (x[[0, 3 + kseg]] + x[[0, 3 + nseg_cd2 + kseg]]) / 2.0,
                    (x[[1, 3 + kseg]] + x[[1, 3 + nseg_cd2 + kseg]]) / 2.0,
                    (x[[2, 3 + kseg]] + x[[2, 3 + nseg_cd2 + kseg]]) / 2.0,
                ];
                let dirx = [
                    x[[0, 3 + nseg_cd2 + kseg]] - midx[0],
                    x[[1, 3 + nseg_cd2 + kseg]] - midx[1],
                    x[[2, 3 + nseg_cd2 + kseg]] - midx[2],
                ];

                let ctype = secstr_cd2[kseg];

                // Grid position
                let gx0 = fung(midx[0]);
                let gy0 = fung(midx[1]);
                let gz0 = fung(midx[2]);

                // Search grid neighborhood
                for entry in grid.lookup(gx0, gy0, gz0, 2) {
                    let cur_a = entry.aseg;
                    let cur_b = entry.bseg;
                    let cur_c = entry.cseg;

                    // Filter: secondary structure type match
                    if secstr_cd1[cur_a] != atype {
                        continue;
                    }
                    if secstr_cd1[cur_c] != ctype {
                        continue;
                    }

                    // Topology filters
                    if less {
                        if cur_a >= cur_c {
                            continue;
                        }
                    } else if cur_a <= cur_c {
                        continue;
                    }

                    if lest {
                        if cur_a >= cur_b {
                            continue;
                        }
                    } else if cur_a <= cur_b {
                        continue;
                    }

                    // Midpoint distance filter (< 4.0 Å)
                    let mid_entry = [
                        (entry.link_from[0] + entry.link_to[0]) / 2.0,
                        (entry.link_from[1] + entry.link_to[1]) / 2.0,
                        (entry.link_from[2] + entry.link_to[2]) / 2.0,
                    ];
                    let dir_entry = [
                        entry.link_to[0] - mid_entry[0],
                        entry.link_to[1] - mid_entry[1],
                        entry.link_to[2] - mid_entry[2],
                    ];

                    let dist = distance_3d(&mid_entry, &midx);
                    if dist > 4.0 {
                        continue;
                    }

                    // Direction cosine filter (> 0.5)
                    if cosi(&dir_entry, &dirx) < 0.5 {
                        continue;
                    }

                    // Vote
                    *count.entry((cur_a, cur_b)).or_insert(0) += 1;
                }
            }

            // Find best (useg, vseg) for this cd2 pair
            // Fortran iterates over protein_nseg bounds (cd2's nseg)
            let protein_nseg = nseg_cd1;
            let mut local_best = 0usize;
            let mut local_useg = 0usize;
            let mut local_vseg = 0usize;
            for lseg in 0..protein_nseg {
                for kseg_inner in 0..protein_nseg {
                    let c = count.get(&(lseg, kseg_inner)).copied().unwrap_or(0);
                    if c > local_best {
                        local_best = c;
                        local_useg = lseg;
                        local_vseg = kseg_inner;
                    }
                }
            }

            if local_best > protcount {
                protcount = local_best;
                best = Some(CompareResult {
                    iseg,
                    jseg,
                    useg: local_useg,
                    vseg: local_vseg,
                    count: local_best,
                });
            }
        }
    }

    best
}

#[inline]
fn distance_3d(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

#[inline]
fn cosi(z: &[f64; 3], x: &[f64; 3]) -> f64 {
    let lz = (z[0] * z[0] + z[1] * z[1] + z[2] * z[2]).sqrt();
    let lx = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
    let norm = (lz * lx).max(1e-6);
    (z[0] * x[0] + z[1] * x[1] + z[2] * x[2]) / norm
}
