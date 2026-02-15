//! Protein I/O for PARSI — reads .dat files including domain tree, distance matrices.

use ndarray::Array2;
use crate::numerics::rounding::nint;
use super::{MAXRES1, NUL};

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
    pub node_child: Vec<i32>,        // flat [child0_0, child0_1, child1_0, child1_1, ...]
    pub domns: Vec<i32>,             // number of segments per domain
    pub domseglist: Vec<i32>,        // flat [seg indices], indexed by domseglist[j * (ndom+1) + idom]
    pub domseglist_stride: usize,    // max_segs_per_dom (row stride)
}

impl ParsiProtein {
    /// Get node_child[side, idom] (side: 0=left, 1=right).
    #[inline]
    pub fn child(&self, side: usize, idom: usize) -> i32 {
        self.node_child[side * (self.ndom + 1) + idom]
    }

    /// Set node_child[side, idom].
    #[inline]
    pub fn set_child(&mut self, side: usize, idom: usize, val: i32) {
        self.node_child[side * (self.ndom + 1) + idom] = val;
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
pub fn parsireadproteindata(filepath: &str) -> Option<ParsiProtein> {
    let content = std::fs::read_to_string(filepath).ok()?;
    let lines: Vec<&str> = content.lines().collect();
    let mut idx = 0;

    if lines.is_empty() {
        return None;
    }

    // Section 1: Header (format 500: 10x,4i5,2x,200a1)
    let header = lines[idx]; idx += 1;
    if header.len() < 32 {
        return None;
    }
    let nres: usize = header[10..15].trim().parse().ok()?;
    let nseg: usize = header[15..20].trim().parse().ok()?;
    let na: usize = header[20..25].trim().parse().ok()?;
    let nb: usize = header[25..30].trim().parse().ok()?;
    let secstr_str = &header[32..32 + nseg.min(header.len() - 32)];
    let secstr: Vec<char> = secstr_str.chars().collect();
    let code = header[5..10].trim().to_string();

    // Section 2: Segment ranges (format 510: 6i10)
    let mut segmentrange = vec![0i32; nseg * 2];
    let mut checkrange = vec![0i32; nseg * 2];
    let mut checkx = vec![0i32; nseg];
    for _iseg in 0..nseg {
        if idx >= lines.len() { return None; }
        let line = lines[idx]; idx += 1;
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
        if idx >= lines.len() { return None; }
        let line = lines[idx]; idx += 1;
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

    // Section 4: Domain tree
    if idx >= lines.len() { return None; }
    let domline = lines[idx]; idx += 1;
    let ndom: usize = domline[10..15].trim().parse().unwrap_or(0);

    let alloc = ndom + 100; // extra space for treehack growth
    let mut node_child = vec![0i32; 2 * alloc];
    let mut domns = vec![0i32; alloc];
    let mut domseglist_lists: Vec<Vec<i32>> = vec![Vec::new(); alloc];

    for _ in 0..ndom {
        if idx >= lines.len() { break; }
        let line = lines[idx]; idx += 1;
        if line.len() < 19 { continue; }
        let idom: usize = line[0..4].trim().parse().unwrap_or(0);
        let _nt = line.chars().nth(5).unwrap_or(' ');
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

        if idom < alloc {
            node_child[0 * alloc + idom] = child1;  // left
            node_child[1 * alloc + idom] = child2;  // right
            domns[idom] = dns;
            domseglist_lists[idom] = segs;
        }
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

    Some(ParsiProtein {
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
