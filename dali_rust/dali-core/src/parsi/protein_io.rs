//! Protein I/O for PARSI — reads .dat files including domain tree, distance matrices.

use ndarray::Array2;
use crate::numerics::rounding::nint;
use super::{MAXDOM, MAXRES1, MAXSEG, NUL};

/// Full protein data for PARSI (includes domain tree).
#[derive(Debug, Clone)]
pub struct ParsiProtein {
    pub code: String,
    pub nres: usize,
    pub nseg: usize,
    pub na: usize,
    pub nb: usize,
    pub secstr: Vec<char>,           // 'H' or 'E' per SSE
    pub segmentrange: Vec<i32>,      // flat [s0, e0, s1, e1, ...] 1-based (compressed)
    pub segmentrange0: Vec<i32>,     // flat [s0, e0, ...] 1-based (original)
    pub ca: Array2<f32>,             // (3, nres) CA coordinates (f32 for Fortran compat)
    pub checkrange: Vec<i32>,        // flat [cs0, ce0, cs1, ce1, ...]
    pub checkx: Vec<i32>,
    // Domain tree
    pub ndom: usize,
    pub node_child: Vec<i32>,        // flat [side * ndom_alloc + idom]; right children start at ndom_alloc
    pub domns: Vec<i32>,             // number of segments per domain
    pub domseglist: Vec<i32>,        // flat [seg indices], indexed by domseglist[j * (ndom+1) + idom]
    pub domseglist_stride: usize,    // max_segs_per_dom (row stride)
}

impl ParsiProtein {
    /// Get node_child[side, idom] (side: 0=left, 1=right).
    #[inline]
    pub fn child(&self, side: usize, idom: usize) -> i32 {
        self.node_child[side * self.ndom_alloc() + idom]
    }

    /// Set node_child[side, idom].
    #[inline]
    pub fn set_child(&mut self, side: usize, idom: usize, val: i32) {
        let stride = self.ndom_alloc();
        self.node_child[side * stride + idom] = val;
    }

    /// Get domseglist[j, idom].
    #[inline]
    pub fn seglist(&self, j: usize, idom: usize) -> i32 {
        self.domseglist[j * (self.ndom_alloc()) + idom]
    }

    /// Set domseglist[j, idom].
    #[inline]
    pub fn set_seglist(&mut self, j: usize, idom: usize, val: i32) {
        let stride = self.ndom_alloc();
        self.domseglist[j * stride + idom] = val;
    }

    fn ndom_alloc(&self) -> usize {
        // The allocation stride for node_child / domns / domseglist
        self.node_child.len() / 2
    }

    /// Get segmentrange (start, end) for segment iseg (0-based).
    #[inline]
    pub fn seg_range(&self, iseg: usize) -> (i32, i32) {
        (self.segmentrange[iseg * 2], self.segmentrange[iseg * 2 + 1])
    }

    /// Get segmentrange0 (start, end) for segment iseg (0-based).
    #[inline]
    pub fn seg_range0(&self, iseg: usize) -> (i32, i32) {
        (self.segmentrange0[iseg * 2], self.segmentrange0[iseg * 2 + 1])
    }
}

/// Read .dat file with domain tree (full PARSI format).
///
/// Fortran parsireadproteindata reads section 1 with (10x,4i5,2x,200a1) and then the
/// SSE-level domain tree with (i4,1x,a1,1x,3i4,200i4). Headers are tokenized on
/// whitespace instead, so codes longer than 5 characters work. Files written by the
/// pre-2026-10 Rust writer have no SSE tree (section 2 is a bare header); for those the
/// tree is derived from the residue-level tree in section 3, as puutos.f does.
pub fn parsireadproteindata(filepath: &str) -> Result<ParsiProtein, String> {
    let content = std::fs::read_to_string(filepath).map_err(|e| e.to_string())?;
    let lines: Vec<&str> = content.lines().collect();
    let mut idx = 0;

    // Section 1: Header
    let header = *lines.first().ok_or("empty file")?;
    idx += 1;
    let tok: Vec<&str> = header.split_whitespace().collect();
    if tok.len() < 6 || tok[0] != ">>>>" {
        return Err(format!("bad header '{}'", header));
    }
    let field = |i: usize, name: &str| -> Result<usize, String> {
        tok[i].parse().map_err(|_| format!("bad {} '{}' in header", name, tok[i]))
    };
    let code = tok[1].to_string();
    let nres = field(2, "nres")?;
    let nseg = field(3, "nseg")?;
    let na = field(4, "na")?;
    let nb = field(5, "nb")?;
    if nseg > MAXSEG {
        return Err(format!("nseg {} exceeds MAXSEG {}", nseg, MAXSEG));
    }
    let secstr: Vec<char> = tok.get(6).map(|s| s.chars().take(nseg).collect()).unwrap_or_default();
    if secstr.len() < nseg {
        return Err(format!("secstr has {} chars for {} segments", secstr.len(), nseg));
    }

    // Section 2: Segment ranges (format 510: 6i10)
    let mut segmentrange = vec![0i32; nseg * 2];
    let mut checkrange = vec![0i32; nseg * 2];
    let mut checkx = vec![0i32; nseg];
    for _iseg in 0..nseg {
        let line = *lines.get(idx).ok_or("truncated in segments")?;
        idx += 1;
        let mut vals = Vec::new();
        for k in 0..6 {
            let s = k * 10;
            let e = (s + 10).min(line.len());
            if s >= line.len() { break; }
            if let Ok(v) = line[s..e].trim().parse::<i32>() {
                vals.push(v);
            }
        }
        if vals.len() >= 6 {
            let seg_idx = (vals[0] - 1) as usize; // 1-based to 0-based
            segmentrange[seg_idx * 2] = vals[1];
            segmentrange[seg_idx * 2 + 1] = vals[2];
            checkrange[seg_idx * 2] = vals[3];
            checkrange[seg_idx * 2 + 1] = vals[4];
            checkx[seg_idx] = vals[5];
        }
    }

    // Section 3: CA coordinates (format 520: 10f8.1)
    let total_vals = 3 * nres;
    let mut ca_vals = Vec::with_capacity(total_vals);
    while ca_vals.len() < total_vals {
        let line = *lines.get(idx).ok_or("truncated in CA coordinates")?;
        idx += 1;
        let nfloats = (total_vals - ca_vals.len()).min(10);
        for k in 0..nfloats {
            let s = k * 8;
            let e = (s + 8).min(line.len());
            if s >= line.len() { break; }
            if let Ok(v) = line[s..e].trim().parse::<f32>() {
                ca_vals.push(v);
            }
        }
    }
    // Reshape: input is (nres, 3), we want (3, nres)
    let mut ca = Array2::<f32>::zeros((3, nres));
    for i in 0..nres {
        for d in 0..3 {
            ca[[d, i]] = ca_vals[i * 3 + d];
        }
    }

    // Section 4: SSE-level domain tree, as (child1, child2, seglist) per node
    let domline = *lines.get(idx).ok_or("no domain tree section")?;
    idx += 1;
    let dtok: Vec<&str> = domline.split_whitespace().collect();
    let nodes: Vec<(i32, i32, Vec<i32>)> = if dtok.len() == 3 && dtok[0] == ">>>>" {
        let ndom: usize = dtok[2].parse()
            .map_err(|_| format!("bad ndom in '{}'", domline))?;
        let mut nodes = vec![(0, 0, Vec::new()); ndom];
        for _ in 0..ndom {
            let line = *lines.get(idx).ok_or("truncated in domain tree")?;
            idx += 1;
            if line.len() < 19 { continue; }
            let idom: usize = line[0..4].trim().parse().unwrap_or(0);
            let child1: i32 = line[7..11].trim().parse().unwrap_or(0);
            let child2: i32 = line[11..15].trim().parse().unwrap_or(0);
            let dns: i32 = line[15..19].trim().parse().unwrap_or(0);

            let mut segs = Vec::new();
            let mut pos = 19;
            for _ in 0..dns {
                let end = (pos + 4).min(line.len());
                if pos >= line.len() { break; }
                if let Ok(v) = line[pos..end].trim().parse::<i32>() {
                    segs.push(v);
                }
                pos += 4;
            }
            if idom >= 1 && idom <= ndom {
                nodes[idom - 1] = (child1, child2, segs);
            }
        }
        nodes
    } else {
        // No SSE tree in the file: derive it from the residue-level tree.
        let prot = crate::io::read_dat(filepath).map_err(|e| e.to_string())?;
        crate::io::derive_sse_tree(prot.nres, &prot.segments, &prot.domain_tree)
            .into_iter()
            .map(|n| (n.left_child as i32, n.right_child as i32,
                      n.sses.iter().map(|&s| s as i32).collect()))
            .collect()
    };

    let ndom = nodes.len();
    if ndom > MAXDOM {
        return Err(format!("ndom {} exceeds MAXDOM {}", ndom, MAXDOM));
    }

    // treehack adds at most 2 nodes per SSE beyond the first in each leaf
    let alloc = ndom + 2 * nseg + 100;
    let mut node_child = vec![0i32; 2 * alloc];
    let mut domns = vec![0i32; alloc];
    let mut domseglist_lists: Vec<Vec<i32>> = vec![Vec::new(); alloc];
    for (i, (child1, child2, segs)) in nodes.into_iter().enumerate() {
        let idom = i + 1;
        node_child[0 * alloc + idom] = child1;  // left
        node_child[1 * alloc + idom] = child2;  // right
        domns[idom] = segs.len() as i32;
        domseglist_lists[idom] = segs;
    }

    // Convert domseglist to flat array
    let max_segs = domseglist_lists.iter()
        .map(|sl| sl.len())
        .max()
        .unwrap_or(0)
        .max(nseg);
    let mut domseglist = vec![0i32; max_segs * alloc];
    for idom in 0..alloc {
        for (j, &seg) in domseglist_lists[idom].iter().enumerate() {
            domseglist[j * alloc + idom] = seg;
        }
    }

    Ok(ParsiProtein {
        code, nres, nseg, na, nb, secstr,
        segmentrange: segmentrange.clone(),
        segmentrange0: segmentrange,
        ca, checkrange, checkx,
        ndom,
        node_child,
        domns,
        domseglist,
        domseglist_stride: max_segs,
    })
}

/// Split multi-segment leaf nodes into binary tree.
pub fn treehack(prot: &mut ParsiProtein) {
    let mut ndom = prot.ndom;
    let alloc = prot.ndom_alloc();

    let mut idom = 0;
    while idom < ndom {
        idom += 1; // 1-based
        if prot.domns[idom] > 1 && prot.child(0, idom) == 0 {
            // Leaf with multiple segments — split
            if ndom + 2 >= alloc {
                // Need to grow arrays — for now just skip
                break;
            }

            // Left child: first segment
            ndom += 1;
            prot.set_child(0, idom, ndom as i32);
            prot.domns[ndom] = 1;
            let first_seg = prot.seglist(0, idom);
            prot.set_seglist(0, ndom, first_seg);
            prot.set_child(0, ndom, 0);
            prot.set_child(1, ndom, 0);

            // Right child: remaining segments
            ndom += 1;
            prot.set_child(1, idom, ndom as i32);
            prot.domns[ndom] = prot.domns[idom] - 1;
            for j in 1..prot.domns[idom] as usize {
                let seg = prot.seglist(j, idom);
                prot.set_seglist(j - 1, ndom, seg);
            }
            prot.set_child(0, ndom, 0);
            prot.set_child(1, ndom, 0);
        }
    }
    prot.ndom = ndom;
}

/// Remap CA coordinates to packed indices (remove inter-segment gaps).
pub fn compressca(prot: &mut ParsiProtein) {
    let nseg = prot.nseg;
    let segmentrange0 = &prot.segmentrange0;

    // Build compressed ranges
    let mut sr = vec![0i32; nseg * 2];
    let mut l = 0i32;
    for iseg in 0..nseg {
        l += 1;
        sr[iseg * 2] = l;
        l += segmentrange0[iseg * 2 + 1] - segmentrange0[iseg * 2];
        sr[iseg * 2 + 1] = l;
    }

    let new_nres = if nseg > 0 { sr[(nseg - 1) * 2 + 1] as usize } else { 0 };
    if new_nres > MAXRES1 || new_nres == 0 {
        prot.nres = 0;
        return;
    }

    // Copy coordinates in compressed order
    let mut cx = Array2::<f32>::zeros((3, new_nres));
    for iseg in 0..nseg {
        let cs = sr[iseg * 2];
        let ce = sr[iseg * 2 + 1];
        let os = segmentrange0[iseg * 2];
        for ires in cs..=ce {
            let a = (os + ires - cs) as usize; // 1-based original
            let ci = (ires - 1) as usize;       // 0-based compressed
            if a >= 1 && a <= prot.ca.ncols() && ci < new_nres {
                for d in 0..3 {
                    cx[[d, ci]] = prot.ca[[d, a - 1]];
                }
            }
        }
    }

    prot.nres = new_nres;
    prot.ca = cx;
    prot.segmentrange = sr;
}

/// Compute CA distance matrix for protein 1.
/// dist[i*nres+j] = nint(10*distance), capped at 1000, diagonal=1.
/// Uses f32 arithmetic to match Fortran.
pub fn getdist(prot: &ParsiProtein) -> Vec<i16> {
    let nres = prot.nres;
    let ca = &prot.ca;
    let mut dist = vec![0i16; nres * nres];

    for i in 0..nres {
        dist[i * nres + i] = 1;
        for j in (i + 1)..nres {
            let dx = ca[[0, i]] - ca[[0, j]];
            let dy = ca[[1, i]] - ca[[1, j]];
            let dz = ca[[2, i]] - ca[[2, j]];
            let d = (dx * dx + dy * dy + dz * dz).sqrt();
            let v = nint((10.0f32 * d) as f64) as i16;
            let v = v.min(1000);
            dist[i * nres + j] = v;
            dist[j * nres + i] = v;
        }
    }

    dist
}

/// Cap distances at 400 (= 40.0 Å) in-place.
pub fn hackdist(dist: &mut Vec<i16>, nres: usize) {
    for i in 0..nres * nres {
        if dist[i] > 400 {
            dist[i] = 400;
        }
    }
}

/// Compute binned distance matrix for protein 2.
pub fn getdist2(ca: &Array2<f32>, nres: usize) -> Vec<i16> {
    let mut dist2 = vec![160i16; nres * nres];

    for i in 0..nres {
        dist2[i * nres + i] = 1;
        for j in (i + 1)..nres {
            let dx = ca[[0, i]] - ca[[0, j]];
            let dy = ca[[1, i]] - ca[[1, j]];
            let dz = ca[[2, i]] - ca[[2, j]];
            let x = (dx * dx + dy * dy + dz * dz).sqrt();
            let v = if x <= 10.0f32 {
                nint((x / 0.1f32) as f64) as i16
            } else if x <= 20.0f32 {
                100 + nint(((x - 10.0f32) / 0.4f32) as f64) as i16
            } else if x < 55.0f32 {
                125 + nint((x - 20.0f32) as f64) as i16
            } else {
                160
            };
            dist2[i * nres + j] = v;
            dist2[j * nres + i] = v;
        }
    }

    dist2
}

/// Compute cumulative distance sums.
/// dist2sum[i*(nres+1)+j] = sum of max(100, nint(1000*dist(ca[i], ca[k]))) for k=0..j-1.
pub fn getdist2sum(ca: &Array2<f32>, nres: usize) -> Vec<i32> {
    let stride = nres + 1;
    let mut dist2sum = vec![0i32; nres * stride];

    for i in 0..nres {
        for j in 0..nres {
            let dx = ca[[0, i]] - ca[[0, j]];
            let dy = ca[[1, i]] - ca[[1, j]];
            let dz = ca[[2, i]] - ca[[2, j]];
            let d = (dx * dx + dy * dy + dz * dz).sqrt();
            let v = nint((1000.0f32 * d) as f64).max(100) as i32;
            dist2sum[i * stride + j + 1] = dist2sum[i * stride + j] + v;
        }
    }

    dist2sum
}
