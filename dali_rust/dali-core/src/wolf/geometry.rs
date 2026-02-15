use ndarray::Array2;

/// Compute SSE midpoints and direction vectors from CA coordinates.
///
/// For each SSE, split residues into N-half and C-half:
///   mid = average of (N-half center + C-half center) / 2
///   dir = C-half center - N-half center
///
/// Returns (midpoint, direction) both (3, nseg).
pub fn vec_sse(
    ca: &Array2<f64>,
    _nres: usize,
    nseg: usize,
    seg_starts: &[u32],
    seg_ends: &[u32],
) -> (Array2<f64>, Array2<f64>) {
    let mut midpoint = Array2::zeros((3, nseg));
    let mut direction = Array2::zeros((3, nseg));

    for iseg in 0..nseg {
        let left = seg_starts[iseg] as usize; // 1-based
        let rite = seg_ends[iseg] as usize; // 1-based
        let mid_res = (left + rite) / 2; // Fortran integer division
        let l = (mid_res - left + 1) as f64; // N-half count
        let r = (rite - mid_res + 1) as f64; // C-half count (mid_res shared)

        // N-half: residues left..mid_res (inclusive, 1-based)
        let mut nmid = [0.0f64; 3];
        for ires in left..=mid_res {
            for k in 0..3 {
                nmid[k] += ca[[k, ires - 1]] / l;
            }
        }

        // C-half: residues mid_res..rite (inclusive, 1-based)
        let mut cmid = [0.0f64; 3];
        for ires in mid_res..=rite {
            for k in 0..3 {
                cmid[k] += ca[[k, ires - 1]] / r;
            }
        }

        for k in 0..3 {
            midpoint[[k, iseg]] = (nmid[k] + cmid[k]) / 2.0;
            direction[[k, iseg]] = cmid[k] - nmid[k];
        }
    }

    (midpoint, direction)
}

/// Transform coordinates to canonical frame.
///
/// Puts x[:,0] at origin, x[:,1] along +y axis, x[:,2] in positive yz-plane.
/// Modifies x in-place.
pub fn twist(x: &mut Array2<f64>, looplen: usize) {
    let pi = std::f64::consts::PI;

    // Set origin: subtract x[:,0] from all points
    let origin = [x[[0, 0]], x[[1, 0]], x[[2, 0]]];
    for i in (0..looplen).rev() {
        x[[0, i]] -= origin[0];
        x[[1, i]] -= origin[1];
        x[[2, i]] -= origin[2];
    }

    // Rotate around x-axis to put x[:,1] in xy-plane
    // If y ~= 0, first rotate 90° around z
    if x[[1, 1]].abs() < 1e-6 {
        for i in 1..looplen {
            let y0 = -x[[1, i]];
            let y1 = x[[0, i]];
            let y2 = x[[2, i]];
            x[[0, i]] = y0;
            x[[1, i]] = y1;
            x[[2, i]] = y2;
        }
    }

    // If y ~= 0 and x ~= 0, rotate 90° around x
    if x[[1, 1]].abs() < 1e-6 && x[[0, 1]].abs() < 1e-6 {
        for i in 1..looplen {
            let y1 = -x[[2, i]];
            let y2 = x[[1, i]];
            let y0 = x[[0, i]];
            x[[0, i]] = y0;
            x[[1, i]] = y1;
            x[[2, i]] = y2;
        }
    }

    // Rotate around x-axis
    let u = if x[[1, 1]].abs() > 1e-6 {
        (x[[2, 1]] / x[[1, 1]]).atan()
    } else {
        0.0
    };
    let sinu = u.sin();
    let cosu = u.cos();
    for i in 1..looplen {
        let y1 = cosu * x[[1, i]] + sinu * x[[2, i]];
        let y2 = -sinu * x[[1, i]] + cosu * x[[2, i]];
        x[[1, i]] = y1;
        x[[2, i]] = y2;
    }

    // Rotate around z-axis
    let u = if x[[1, 1]].abs() > 1e-6 {
        -(x[[0, 1]] / x[[1, 1]]).atan()
    } else {
        0.0
    };
    let u = if x[[1, 1]] < 0.0 { pi + u } else { u };
    let sinu = u.sin();
    let cosu = u.cos();
    for i in 1..looplen {
        let y0 = cosu * x[[0, i]] + sinu * x[[1, i]];
        let y1 = -sinu * x[[0, i]] + cosu * x[[1, i]];
        x[[0, i]] = y0;
        x[[1, i]] = y1;
    }

    // Rotate around y-axis to put x[:,2] in yz-plane
    // If z ~= 0, first rotate 90° around y
    if x[[2, 2]].abs() < 1e-6 {
        for i in 1..looplen {
            let y0 = x[[2, i]];
            let y2 = -x[[0, i]];
            x[[0, i]] = y0;
            x[[2, i]] = y2;
        }
    }

    let u = if x[[2, 2]] != 0.0 {
        (x[[0, 2]] / x[[2, 2]]).atan()
    } else {
        0.0
    };
    let sinu = u.sin();
    let cosu = u.cos();
    for i in 2..looplen {
        let y2 = cosu * x[[2, i]] + sinu * x[[0, i]];
        let y0 = -sinu * x[[2, i]] + cosu * x[[0, i]];
        x[[0, i]] = y0;
        x[[2, i]] = y2;
    }

    // If z < 0, rotate 180° around y
    if x[[2, 2]] < 0.0 {
        for i in 1..looplen {
            x[[0, i]] = -x[[0, i]];
            x[[2, i]] = -x[[2, i]];
        }
    }
}

/// Build reference frame array for twist().
///
/// Layout:
///   x[:,0] = midpoint[:,iseg]
///   x[:,1] = midpoint[:,iseg] + direction[:,iseg]
///   x[:,2] = midpoint[:,jseg]
///   x[:,3:3+nseg] = midpoint[:,k] - direction[:,k]  for k=0..nseg-1
///   x[:,3+nseg:3+2*nseg] = midpoint[:,k] + direction[:,k]  for k=0..nseg-1
pub fn preparex(
    iseg: usize,
    jseg: usize,
    nseg: usize,
    midpoint: &Array2<f64>,
    direction: &Array2<f64>,
) -> Array2<f64> {
    let n = 3 + 2 * nseg;
    let mut x = Array2::zeros((3, n));

    for k in 0..3 {
        x[[k, 0]] = midpoint[[k, iseg]];
        x[[k, 1]] = midpoint[[k, iseg]] + direction[[k, iseg]];
        x[[k, 2]] = midpoint[[k, jseg]];
    }

    for kseg in 0..nseg {
        for k in 0..3 {
            x[[k, 3 + kseg]] = midpoint[[k, kseg]] - direction[[k, kseg]];
            x[[k, 3 + nseg + kseg]] = midpoint[[k, kseg]] + direction[[k, kseg]];
        }
    }

    x
}

/// Compute pairwise Euclidean distances between SSE midpoints.
pub fn compute_neidist(nseg: usize, midpoint: &Array2<f64>) -> Array2<f64> {
    let mut neidist = Array2::zeros((nseg, nseg));
    for i in 0..nseg {
        for j in 0..nseg {
            if i != j {
                let dx = midpoint[[0, i]] - midpoint[[0, j]];
                let dy = midpoint[[1, i]] - midpoint[[1, j]];
                let dz = midpoint[[2, i]] - midpoint[[2, j]];
                neidist[[i, j]] = (dx * dx + dy * dy + dz * dz).sqrt();
            }
        }
    }
    neidist
}
