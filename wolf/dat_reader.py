"""
Parse DaliLite .dat files into Protein dataclass.

Format (Fortran fixed-width):
  Header: (10x, 4i5, 2x, 200a1)  -> skip 10, 4 ints (5-wide), skip 2, secstr chars
  Segments: (6i10)                -> index, start, end, check_start, check_end, checkx
  CA coords: (10f8.1)            -> 10 floats per line, 8 chars wide
  Domain info: follows CA coords  -> read but not used
"""

import numpy as np
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Protein:
    code: str
    nres: int
    nseg: int
    na: int  # number of helices
    nb: int  # number of strands
    secstr: list  # ['H','H','E',...] one per SSE
    segmentrange: np.ndarray  # (2, nseg) — [0,:] = start, [1,:] = end (1-based)
    ca: np.ndarray  # (3, nres) — CA coordinates (0-indexed columns)


def read_dat(filepath):
    """Parse a .dat file into a Protein object."""
    filepath = Path(filepath)
    with open(filepath) as f:
        lines = f.readlines()

    idx = 0

    # Header line: format (10x, 4i5, 2x, 200a1)
    header = lines[idx]
    idx += 1
    # Skip 10 chars, then read 4 integers (5 chars each)
    nres = int(header[10:15])
    nseg = int(header[15:20])
    na = int(header[20:25])
    nb = int(header[25:30])
    # Skip 2 chars (positions 30-31), then read secstr chars
    secstr_str = header[32:32 + nseg]
    secstr = list(secstr_str)

    # Extract code from the header (chars 5-9, typically)
    code = header[5:10].strip()

    # Segment ranges: format (6i10) — one line per segment
    segmentrange = np.zeros((2, nseg), dtype=np.int32)
    for iseg in range(nseg):
        line = lines[idx]
        idx += 1
        # 6 integers, 10 chars each: index, start, end, check_start, check_end, checkx
        vals = []
        for k in range(6):
            s = line[k * 10:(k + 1) * 10]
            vals.append(int(s))
        segmentrange[0, iseg] = vals[1]  # start (1-based)
        segmentrange[1, iseg] = vals[2]  # end (1-based)

    # CA coordinates: format (10f8.1)
    # Total values: 3 * nres (x,y,z for each residue)
    total_vals = 3 * nres
    ca_vals = []
    while len(ca_vals) < total_vals:
        line = lines[idx]
        idx += 1
        # 10 floats per line, 8 chars each
        nfloats = min(10, total_vals - len(ca_vals))
        for k in range(nfloats):
            s = line[k * 8:(k + 1) * 8]
            ca_vals.append(float(s))

    # Reshape: ca(j, i) in Fortran = ca[j-1, i-1] in Python
    # Layout in file: ca(1,1), ca(2,1), ca(3,1), ca(1,2), ca(2,2), ca(3,2), ...
    ca = np.array(ca_vals, dtype=np.float64).reshape(nres, 3).T  # (3, nres)

    return Protein(
        code=code,
        nres=nres,
        nseg=nseg,
        na=na,
        nb=nb,
        secstr=secstr,
        segmentrange=segmentrange,
        ca=ca,
    )
