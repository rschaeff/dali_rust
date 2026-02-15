"""
DALICON output formatting: convert ali1 to WOLFITZ lines.

Translated from gagadccp.f output(), getali0(), getblocks(), filter().
"""

import numpy as np


def getali0(ali1, blocksize, nres1):
    """
    Uncompress alignment from blocksize-compressed to per-residue.

    With blocksize=1 (as called from dowork_dalicon), this just copies ali1.

    Translated from gagadccp.f getali0().

    Args:
        ali1: (nres1+1,) int16 alignment array (1-based)
        blocksize: block size (1 for DALICON output)
        nres1: number of residues

    Returns:
        a1: (nres1+1,) int alignment array (1-based)
    """
    a1 = np.zeros(nres1 + 1, dtype=np.int32)
    for j in range(1, nres1 + 1):
        i1 = int(ali1[j])
        if i1 != 0:
            for l in range(blocksize):
                if j + l <= nres1:
                    a1[j + l] = i1 + l
    return a1


def getblocks(ali1, blocksize, nres1):
    """
    Extract block assignments from alignment.

    Translated from dimple.f getblocks() — blocksize=1 special case.

    Args:
        ali1: (nres1+1,) int16 alignment (1-based)
        blocksize: block size (1 for DALICON output)
        nres1: number of residues

    Returns:
        block: (nres1+1,) int block assignments (1-based)
        nblock: number of blocks
    """
    block = np.zeros(nres1 + 1, dtype=np.int32)
    nblock = 0

    if blocksize > 1:
        for i in range(1, nres1 - blocksize + 2):
            if ali1[i] != 0:
                if block[i] == 0:
                    nblock += 1
                for j in range(blocksize):
                    block[i + j] = nblock
    else:
        # blocksize=1 special case
        if ali1[1] != 0:
            nblock = 1
            block[1] = 1
        for i in range(2, nres1 + 1):
            if ali1[i] != 0:
                # Two SEPARATE if checks (not elif) — matches Fortran exactly
                # Special case: N-term of second sequence
                if ali1[i] == 1:
                    nblock += 1
                if (ali1[i] != ali1[i - 1] + 1) and (ali1[i] != ali1[i - 1] - 1):
                    nblock += 1
                block[i] = nblock

    return block, nblock


def filter_score(score, alilen, nblock, nres1, nres2, minblocks=4):
    """
    Filter: set score to 0 if too few blocks AND too few residues.

    Translated from gagadccp.f filter().

    Args:
        score: alignment score
        alilen: alignment length (residues)
        nblock: number of blocks
        nres1, nres2: sequence lengths
        minblocks: minimum block count

    Returns:
        filtered score (0.0 if filtered out)
    """
    minlen = min(nres1 // 2, nres2 // 2, 30)
    if alilen < minlen and nblock < minblocks:
        return 0.0
    return score


def format_wolfitz(cd1, cd2, nblock, ib1, ib2, a1):
    """
    Format a WOLFITZ output line.

    The Fortran format is:
        write(*,*) 'WOLFITZ ',code1,chainid1,code2,chainid2,
            nblock,(ib1(j),ib2(j),j=1,nblock),(a1(ib1(j)),a1(ib2(j)),j=1,nblock)

    Args:
        cd1, cd2: structure codes (e.g. "101mA")
        nblock: number of blocks
        ib1: block start positions (1-based, list)
        ib2: block end positions (1-based, list)
        a1: per-residue alignment array (1-based)

    Returns:
        WOLFITZ line string
    """
    code1 = cd1[:4]
    chainid1 = cd1[4] if len(cd1) > 4 else ' '
    code2 = cd2[:4]
    chainid2 = cd2[4] if len(cd2) > 4 else ' '

    # Fortran list-directed output with character*4 and character*1
    parts = [f' WOLFITZ {code1}{chainid1}{code2}{chainid2}']
    vals = [nblock]
    for j in range(nblock):
        vals.append(ib1[j])
        vals.append(ib2[j])
    for j in range(nblock):
        vals.append(int(a1[ib1[j]]))
        vals.append(int(a1[ib2[j]]))
    parts.extend(f'{v:>12d}' for v in vals)
    return ''.join(parts)


def produce_output(cd1, cd2, nres1, nres2, ali1, score_val):
    """
    Produce WOLFITZ output line(s) from alignment.

    Translated from gagadccp.f output() (simplified for popsize=1).

    Args:
        cd1, cd2: structure codes
        nres1, nres2: sequence lengths
        ali1: (nres1+1,) int16 alignment (1-based)
        score_val: alignment score

    Returns:
        list of WOLFITZ line strings (0 or 1 lines)
    """
    blocksize = 1  # DALICON uses blocksize=1 for output

    # getali0: uncompress
    a1 = getali0(ali1, blocksize, nres1)

    # getblocks
    block, nblock = getblocks(ali1, blocksize, nres1)

    # Compute alignment length
    lenali = 0
    for j in range(1, nres1 + 1):
        if a1[j] != 0:
            lenali += 1

    # Filter
    filtered_score = filter_score(score_val, lenali, nblock, nres1, nres2)
    if filtered_score <= 0.0:
        return []

    # Build block boundaries (ib1[j] = first residue of block j, ib2[j] = last)
    ib1 = [0] * (nblock + 1)
    ib2 = [0] * (nblock + 1)
    for j in range(1, nres1 + 1):
        if block[j] > 0:
            b = block[j]
            if ib1[b] == 0:
                ib1[b] = j
            ib2[b] = j

    # Convert to 0-based lists for format_wolfitz
    ib1_list = [ib1[b] for b in range(1, nblock + 1)]
    ib2_list = [ib2[b] for b in range(1, nblock + 1)]

    line = format_wolfitz(cd1, cd2, nblock, ib1_list, ib2_list, a1)
    return [line]
