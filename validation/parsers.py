"""
Parsers for DaliLite intermediate file formats.

Each parser reads a specific format into Python data structures
that can be compared with tolerance.
"""

import re
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class WolfitzHit:
    """A single WOLFITZ line: initial SSE-based alignment."""
    cd1: str
    cd2: str
    nblock: int
    # Residue ranges: list of (start, end) for each block in cd1 and cd2
    blocks_cd1: list  # [(start, end), ...]
    blocks_cd2: list  # [(start, end), ...]
    raw_values: list  # All integer values after nblock


@dataclass
class DCCPBlock:
    """A single DCCP alignment block."""
    cd1: str
    cd2: str
    score: float
    zscore: float
    rmsd: float
    nblock: int
    # Residue ranges for each block
    ranges_cd1: list  # [(start, end), ...]
    ranges_cd2: list  # [(start, end), ...]


@dataclass
class RefineLine:
    """A 'refine' output line from PARSI."""
    cd1cd2: str  # 10-char string
    idom: int
    score: int
    nseg: int
    # Segment ranges: (a1, a2) for cd1, (b1, b2) for cd2
    ranges_cd1: list
    ranges_cd2: list


@dataclass
class Filter95Line:
    """A FILTER95 output line: filtered/refined alignment."""
    cd1cd2: str  # 10-char string
    cd1: str
    cd2: str
    idom: int
    score: int
    zscore: float
    nres: int
    # Per-residue alignment: (a1, a2) pairs for cd1 and cd2
    ranges_cd1: list  # [(a1, a2), ...]
    ranges_cd2: list  # [(b1, b2), ...]


@dataclass
class DaliconInput:
    """A DALICON input record: pre-alignment to refine."""
    cd1: str
    cd2: str
    nblock: int
    ranges_cd1: list  # [(start, end), ...]
    ranges_cd2: list  # [(start, end), ...]


def parse_wolfitz(filepath):
    """Parse WOLFITZ output file into list of WolfitzHit."""
    hits = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line.startswith('WOLFITZ'):
                continue
            parts = line.split()
            # WOLFITZ cd1cd2 nblock values...
            cd1cd2 = parts[1]
            cd1 = cd1cd2[:5].strip()
            cd2 = cd1cd2[5:].strip()
            nblock = int(parts[2])
            values = [int(x) for x in parts[3:]]
            hits.append(WolfitzHit(
                cd1=cd1, cd2=cd2, nblock=nblock,
                blocks_cd1=[], blocks_cd2=[],
                raw_values=values,
            ))
    return hits


def parse_dccp(filepath):
    """Parse DCCP file into list of DCCPBlock."""
    blocks = []
    with open(filepath) as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith('DCCP'):
            # DCCP   1   score zscore   0   rmsd   0   nblock   cd1 cd2
            parts = line.split()
            score = float(parts[2])
            zscore = float(parts[3])
            rmsd = float(parts[5])
            nblock = int(parts[7])
            cd1 = parts[8]
            cd2 = parts[9]

            # Next line should be "alignment"
            i += 1
            if i < len(lines) and 'alignment' in lines[i]:
                i += 1

            # Read cd1 ranges (may span multiple lines)
            range_values = []
            while i < len(lines) and not lines[i].strip().startswith('DCCP'):
                vals = lines[i].split()
                if vals:
                    range_values.extend([int(v) for v in vals])
                i += 1
                # After cd1 ranges, next set is cd2 ranges
                # Total values = nblock * 2 (start, end) * 2 (cd1, cd2)
                if len(range_values) >= nblock * 4:
                    break

            # Split into cd1 and cd2 ranges
            ranges_cd1 = []
            ranges_cd2 = []
            if len(range_values) >= nblock * 4:
                for b in range(nblock):
                    ranges_cd1.append((range_values[b * 2], range_values[b * 2 + 1]))
                offset = nblock * 2
                for b in range(nblock):
                    ranges_cd2.append((range_values[offset + b * 2],
                                       range_values[offset + b * 2 + 1]))

            blocks.append(DCCPBlock(
                cd1=cd1, cd2=cd2, score=score, zscore=zscore,
                rmsd=rmsd, nblock=nblock,
                ranges_cd1=ranges_cd1, ranges_cd2=ranges_cd2,
            ))
            continue
        i += 1

    return blocks


def parse_refine_lines(filepath):
    """Parse FILTER95/PARSI refine output."""
    results = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line.startswith('refine'):
                continue
            # refine<cd1cd2> idom score nseg a1 a2 ... b1 b2 ...
            parts = line.split()
            cd1cd2 = parts[0][6:]  # strip 'refine'
            idom = int(parts[1])
            score = int(parts[2])
            nseg = int(parts[3])
            values = [int(x) for x in parts[4:]]

            ranges_cd1 = []
            ranges_cd2 = []
            for s in range(nseg):
                ranges_cd1.append((values[s * 2], values[s * 2 + 1]))
            offset = nseg * 2
            for s in range(nseg):
                ranges_cd2.append((values[offset + s * 2],
                                   values[offset + s * 2 + 1]))

            results.append(RefineLine(
                cd1cd2=cd1cd2, idom=idom, score=score, nseg=nseg,
                ranges_cd1=ranges_cd1, ranges_cd2=ranges_cd2,
            ))
    return results


def parse_filter95(filepath):
    """Parse FILTER95 output: zscore cd1cd2 idom score nres ranges..."""
    results = []
    with open(filepath) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 5:
                continue
            try:
                zscore = float(parts[0])
            except ValueError:
                continue
            cd1cd2 = parts[1]
            cd1 = cd1cd2[:5].strip()
            cd2 = cd1cd2[5:].strip()
            idom = int(parts[2])
            score = int(parts[3])
            nres = int(parts[4])
            values = [int(x) for x in parts[5:]]

            ranges_cd1 = []
            ranges_cd2 = []
            for r in range(nres):
                if r * 2 + 1 < len(values):
                    ranges_cd1.append((values[r * 2], values[r * 2 + 1]))
            offset = nres * 2
            for r in range(nres):
                if offset + r * 2 + 1 < len(values):
                    ranges_cd2.append((values[offset + r * 2],
                                       values[offset + r * 2 + 1]))

            results.append(Filter95Line(
                cd1cd2=cd1cd2, cd1=cd1, cd2=cd2,
                idom=idom, score=score, zscore=zscore, nres=nres,
                ranges_cd1=ranges_cd1, ranges_cd2=ranges_cd2,
            ))
    return results


def parse_dalicon_input(filepath):
    """Parse DALICON input format.

    Format is grouped by cd1 with END delimiters:
        END
        cd1
        cd2*
        nblock
          range values (may span multiple lines)

        cd2*
        nblock
          range values
        ...
        END
        cd1
        ...
    """
    records = []
    with open(filepath) as f:
        lines = [l.rstrip() for l in f.readlines()]

    i = 0
    cd1 = None
    while i < len(lines):
        line = lines[i].strip()

        # END marks boundary between cd1 groups
        if line == 'END':
            i += 1
            # Next non-empty line after END is the new cd1
            while i < len(lines) and not lines[i].strip():
                i += 1
            if i < len(lines) and lines[i].strip() != 'END':
                cd1 = lines[i].strip().rstrip('*')
                i += 1
            continue

        # Skip blank lines
        if not line:
            i += 1
            continue

        # If we have no cd1 context yet, this line is cd1
        if cd1 is None:
            cd1 = line.rstrip('*')
            i += 1
            continue

        # Otherwise, this should be a cd2 name (possibly with * suffix)
        # Check: if the line looks like a name (not all digits)
        if not line.replace('-', '').replace('.', '').isdigit():
            cd2 = line.rstrip('*')
            i += 1
            # Next line should be nblock
            while i < len(lines) and not lines[i].strip():
                i += 1
            if i >= len(lines):
                break
            nblock = int(lines[i].strip())
            i += 1

            # Read range values (nblock*4 integers: nblock pairs for cd1 + nblock pairs for cd2)
            values = []
            while len(values) < nblock * 4 and i < len(lines):
                vline = lines[i].strip()
                if not vline or vline == 'END':
                    break
                try:
                    values.extend([int(x) for x in vline.split()])
                except ValueError:
                    break
                i += 1

            ranges_cd1 = []
            ranges_cd2 = []
            for b in range(nblock):
                if b * 2 + 1 < len(values):
                    ranges_cd1.append((values[b * 2], values[b * 2 + 1]))
            offset = nblock * 2
            for b in range(nblock):
                if offset + b * 2 + 1 < len(values):
                    ranges_cd2.append((values[offset + b * 2],
                                       values[offset + b * 2 + 1]))

            records.append(DaliconInput(
                cd1=cd1, cd2=cd2, nblock=nblock,
                ranges_cd1=ranges_cd1, ranges_cd2=ranges_cd2,
            ))
        else:
            i += 1

    return records


def parse_dat_header(filepath):
    """Parse the header of a .dat file to get structure metadata."""
    with open(filepath) as f:
        line = f.readline()
        # >>>> cd   nres nseg nh ne  secstr...
        parts = line.split()
        cd = parts[1] if len(parts) > 1 else ''
        nres = int(parts[2]) if len(parts) > 2 else 0
        nseg = int(parts[3]) if len(parts) > 3 else 0
        nh = int(parts[4]) if len(parts) > 4 else 0
        ne = int(parts[5]) if len(parts) > 5 else 0
    return {'cd': cd, 'nres': nres, 'nseg': nseg, 'nh': nh, 'ne': ne}
