"""Protein I/O for PARSI — reads .dat files including domain tree section."""

import numpy as np
from dataclasses import dataclass, field
from pathlib import Path
from .parsizes import MAXRES1, NUL


def _nint(x):
    """Fortran-compatible nint: round half away from zero."""
    return int(np.trunc(x + np.copysign(0.5, x)))


@dataclass
class ParsiProtein:
    """Full protein data for PARSI (includes domain tree)."""
    code: str
    nres: int
    nseg: int
    na: int  # number of helices
    nb: int  # number of strands
    secstr: list  # ['H','H','E',...] per SSE
    segmentrange: np.ndarray  # (2, nseg) 1-based, compressed after compressca
    segmentrange0: np.ndarray  # (2, nseg) 1-based, original
    ca: np.ndarray  # (3, nres) CA coordinates
    checkrange: np.ndarray  # (2, nseg) 1-based
    checkx: np.ndarray  # (nseg,)
    # Domain tree
    ndom: int = 0
    node_child: np.ndarray = field(default_factory=lambda: np.zeros((2, 0), dtype=np.int32))
    node_type: list = field(default_factory=list)
    domns: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.int32))
    domseglist: np.ndarray = field(default_factory=lambda: np.zeros((0, 0), dtype=np.int32))


def parsireadproteindata(filepath):
    """Read .dat file with domain tree (full PARSI format).

    Fortran formats:
      500: (10x,4i5,2x,200a1)  — header
      510: (6i10)               — segment ranges
      520: (10f8.1)             — CA coordinates
      500: (10x,4i5,...)        — domain tree header (ndom)
      530: (i4,1x,a1,1x,3i4,200i4) — domain nodes
    """
    filepath = Path(filepath)
    with open(filepath) as f:
        lines = f.readlines()

    idx = 0

    # === Section 1: Header ===
    # Format 500: (10x,4i5,2x,200a1)
    header = lines[idx]; idx += 1
    nres = int(header[10:15])
    nseg = int(header[15:20])
    na = int(header[20:25])
    nb = int(header[25:30])
    secstr_str = header[32:32 + nseg]
    secstr = list(secstr_str)
    code = header[5:10].strip()

    # === Section 2: Segment ranges ===
    # Format 510: (6i10) — index, start, end, check_start, check_end, checkx
    segmentrange = np.zeros((2, nseg), dtype=np.int32)
    checkrange = np.zeros((2, nseg), dtype=np.int32)
    checkx = np.zeros(nseg, dtype=np.int32)
    for iseg in range(nseg):
        line = lines[idx]; idx += 1
        vals = []
        for k in range(6):
            s = line[k * 10:(k + 1) * 10]
            vals.append(int(s))
        seg_idx = vals[0] - 1  # Fortran 1-based index
        segmentrange[0, seg_idx] = vals[1]
        segmentrange[1, seg_idx] = vals[2]
        checkrange[0, seg_idx] = vals[3]
        checkrange[1, seg_idx] = vals[4]
        checkx[seg_idx] = vals[5]

    # === Section 3: CA coordinates ===
    # Format 520: (10f8.1)
    total_vals = 3 * nres
    ca_vals = []
    while len(ca_vals) < total_vals:
        line = lines[idx]; idx += 1
        nfloats = min(10, total_vals - len(ca_vals))
        for k in range(nfloats):
            s = line[k * 8:(k + 1) * 8]
            ca_vals.append(float(s))
    ca = np.array(ca_vals, dtype=np.float32).reshape(nres, 3).T  # (3, nres)

    # === Section 4: Domain tree ===
    # Format 500 again: (10x,4i5,...) but only ndom is read
    domline = lines[idx]; idx += 1
    ndom = int(domline[10:15])

    # Allocate domain arrays (1-indexed in Fortran, we use 0-indexed with +1 size)
    from .parsizes import MAXSEG, MAXDOM
    node_child = np.zeros((2, ndom + 1), dtype=np.int32)  # extra for treehack growth
    node_type = [''] * (ndom + 1)
    domns = np.zeros(ndom + 1, dtype=np.int32)
    # domseglist: (maxseg, ndom) but we'll use lists for flexibility
    domseglist_lists = [[] for _ in range(ndom + 1)]

    # Format 530: (i4,1x,a1,1x,3i4,200i4)
    for _ in range(ndom):
        line = lines[idx]; idx += 1
        idom = int(line[0:4])  # 1-based
        nt = line[5]  # node_type character
        # 3 integers (i4 each) starting at position 7
        child1 = int(line[7:11])
        child2 = int(line[11:15])
        dns = int(line[15:19])
        # Remaining: dns segment indices (i4 each)
        segs = []
        pos = 19
        for j in range(dns):
            segs.append(int(line[pos:pos + 4]))
            pos += 4
        node_child[0, idom] = child1
        node_child[1, idom] = child2
        node_type[idom] = nt
        domns[idom] = dns
        domseglist_lists[idom] = segs

    # Convert domseglist to numpy array
    max_segs_per_dom = max((len(sl) for sl in domseglist_lists), default=0)
    max_segs_per_dom = max(max_segs_per_dom, nseg)  # at least nseg wide
    domseglist = np.zeros((max_segs_per_dom, ndom + 1), dtype=np.int32)
    for idom in range(1, ndom + 1):
        for j, seg in enumerate(domseglist_lists[idom]):
            domseglist[j, idom] = seg

    return ParsiProtein(
        code=code, nres=nres, nseg=nseg, na=na, nb=nb,
        secstr=secstr,
        segmentrange=segmentrange,
        segmentrange0=segmentrange.copy(),
        ca=ca,
        checkrange=checkrange,
        checkx=checkx,
        ndom=ndom,
        node_child=node_child,
        node_type=node_type,
        domns=domns,
        domseglist=domseglist,
    )


def treehack(prot):
    """Split multi-segment leaf nodes into binary tree.

    Modifies prot.ndom, node_child, domns, domseglist in-place.
    Leaf nodes with >1 segment get split: first segment becomes left child,
    remaining segments become right child.
    """
    ndom = prot.ndom
    node_child = prot.node_child
    domns = prot.domns
    domseglist = prot.domseglist

    # May need to grow arrays
    max_alloc = node_child.shape[1]

    idom = 0
    while idom < ndom:
        idom += 1  # 1-based
        if domns[idom] > 1 and node_child[0, idom] == 0:
            # Leaf with multiple segments — split it

            # Ensure space
            if ndom + 2 >= max_alloc:
                new_alloc = max_alloc + 100
                node_child = _grow2d(node_child, (2, new_alloc))
                domns = _grow1d(domns, new_alloc)
                domseglist = _grow2d(domseglist, (domseglist.shape[0], new_alloc))
                prot.node_child = node_child
                prot.domns = domns
                prot.domseglist = domseglist
                max_alloc = new_alloc

            # Left child: first segment
            ndom += 1
            node_child[0, idom] = ndom
            domns[ndom] = 1
            domseglist[0, ndom] = domseglist[0, idom]
            node_child[0, ndom] = 0
            node_child[1, ndom] = 0

            # Right child: remaining segments
            ndom += 1
            node_child[1, idom] = ndom
            domns[ndom] = domns[idom] - 1
            for j in range(1, domns[idom]):
                domseglist[j - 1, ndom] = domseglist[j, idom]
            node_child[0, ndom] = 0
            node_child[1, ndom] = 0

    prot.ndom = ndom


def _grow1d(arr, new_size):
    new = np.zeros(new_size, dtype=arr.dtype)
    new[:len(arr)] = arr
    return new


def _grow2d(arr, new_shape):
    new = np.zeros(new_shape, dtype=arr.dtype)
    new[:arr.shape[0], :arr.shape[1]] = arr
    return new


def compressca(prot):
    """Remap CA coordinates to packed indices (remove inter-segment gaps).

    Sets prot.segmentrange to compressed ranges (1-based).
    prot.segmentrange0 retains original ranges.
    Updates prot.ca and prot.nres.
    """
    nseg = prot.nseg
    segmentrange0 = prot.segmentrange0
    ca = prot.ca

    # Build compressed ranges
    sr = np.zeros((2, nseg), dtype=np.int32)
    l = 0
    for iseg in range(nseg):
        l += 1
        sr[0, iseg] = l
        l += segmentrange0[1, iseg] - segmentrange0[0, iseg]
        sr[1, iseg] = l

    new_nres = sr[1, nseg - 1]
    if new_nres > MAXRES1:
        return False  # overflow

    # Copy coordinates in compressed order
    cx = np.zeros((3, new_nres + 1), dtype=np.float32)  # 1-based indexing
    for iseg in range(nseg):
        for ires in range(sr[0, iseg], sr[1, iseg] + 1):
            a = segmentrange0[0, iseg] + ires - sr[0, iseg]
            cx[:, ires] = ca[:, a - 1]  # a is 1-based, ca is 0-based

    # Overwrite protein data
    prot.nres = new_nres
    # Store as 0-based indexed array: ca[:, 0..nres-1]
    prot.ca = cx[:, 1:new_nres + 1]
    prot.segmentrange = sr


def getdist(prot):
    """Compute CA distance matrix for protein 1.

    Returns int16 array dist[i,j] = nint(10*distance), capped at 1000.
    Uses float32 arithmetic to match Fortran.
    Indices are 0-based (dist[0..nres-1, 0..nres-1]).
    """
    nres = prot.nres
    ca = prot.ca  # (3, nres), 0-based
    dist = np.zeros((nres, nres), dtype=np.int16)

    for i in range(nres):
        dist[i, i] = 1  # minimal distance
        for j in range(i + 1, nres):
            dx = np.float32(ca[0, i]) - np.float32(ca[0, j])
            dy = np.float32(ca[1, i]) - np.float32(ca[1, j])
            dz = np.float32(ca[2, i]) - np.float32(ca[2, j])
            d = np.float32(np.sqrt(np.float32(dx*dx + dy*dy + dz*dz)))
            v = _nint(np.float32(10.0) * d)
            if v > 1000:
                v = 1000
            dist[i, j] = v
            dist[j, i] = v

    return dist


def hackdist(dist):
    """Cap distances at 400 (= 40.0 Å) in-place."""
    mask = dist > 400
    np.clip(dist, None, 400, out=dist)


def getdist2(ca, nres):
    """Compute binned distance matrix for protein 2.

    Binning:
      dist <= 10 Å: nint(dist/0.1)   → bins 1-100
      10 < dist <= 20: 100 + nint((dist-10)/0.4) → bins ~100-125
      20 < dist < 55: 125 + nint(dist-20) → bins ~125-160
      dist >= 55: 160 (default, pre-filled)

    Returns int16 array, 0-based indexing.
    """
    dist2 = np.full((nres, nres), 160, dtype=np.int16)

    for i in range(nres):
        dist2[i, i] = 1
        for j in range(i + 1, nres):
            dx = np.float32(ca[0, i]) - np.float32(ca[0, j])
            dy = np.float32(ca[1, i]) - np.float32(ca[1, j])
            dz = np.float32(ca[2, i]) - np.float32(ca[2, j])
            x = np.float32(np.sqrt(np.float32(dx*dx + dy*dy + dz*dz)))
            if x <= np.float32(10.0):
                dist2[i, j] = _nint(x / np.float32(0.1))
            elif x <= np.float32(20.0):
                dist2[i, j] = 100 + _nint((x - np.float32(10.0)) / np.float32(0.4))
            elif x < np.float32(55.0):
                dist2[i, j] = 125 + _nint(x - np.float32(20.0))
            # else: stays 160
            dist2[j, i] = dist2[i, j]

    return dist2


def getdist2sum(ca, nres):
    """Compute cumulative distance sums.

    dist2sum[i, j] = sum of distanceint4(ca[i], ca[k]) for k=1..j
    distanceint4 = max(100, nint(1000*dist))

    Returns int32 array of shape (nres, nres+1), 0-based.
    dist2sum[i, 0] = 0, dist2sum[i, j] = cumsum through j (1-based j maps to index j).
    """
    dist2sum = np.zeros((nres, nres + 1), dtype=np.int32)

    for i in range(nres):
        for j in range(nres):
            dx = np.float32(ca[0, i]) - np.float32(ca[0, j])
            dy = np.float32(ca[1, i]) - np.float32(ca[1, j])
            dz = np.float32(ca[2, i]) - np.float32(ca[2, j])
            d = np.float32(np.sqrt(np.float32(dx*dx + dy*dy + dz*dz)))
            v = max(100, _nint(np.float32(1000.0) * d))
            dist2sum[i, j + 1] = dist2sum[i, j] + v

    return dist2sum
