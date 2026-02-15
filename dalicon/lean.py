"""
DALICON lean Monte Carlo optimization.

Translated from lean.f, clean.f, and dimple.f.
All arrays use 1-based indexing matching the Fortran original.
"""

import numpy as np
from .scoring import scorefun, getp, _nint


MAXPAIR = 6000


def _index_qsort(d, ve, lo, hi):
    """
    Non-recursive quicksort of index array d by values in ve.
    Sorts d[lo..hi] (1-based) so that ve[d[lo]] <= ve[d[lo+1]] <= ...

    Translated from util.f j_index_Qsort + partition.
    """
    stack = []
    stack.append((lo, hi))

    while stack:
        top, bottom = stack.pop()
        while top < bottom:
            # Partition
            upper = top
            lower = bottom
            save = d[upper]
            while upper != lower:
                while upper < lower and ve[save] <= ve[d[lower]]:
                    lower -= 1
                if upper != lower:
                    d[upper] = d[lower]
                while upper < lower and ve[save] >= ve[d[upper]]:
                    upper += 1
                if upper != lower:
                    d[lower] = d[upper]
            d[upper] = save
            p = upper

            if (p - top) > (bottom - p):
                stack.append((top, p - 1))
                top = p + 1
            else:
                stack.append((p + 1, bottom))
                bottom = p - 1


def initali1(ali, nres):
    """Zero alignment array (1-based, indices 1..nres)."""
    for i in range(1, nres + 1):
        ali[i] = 0


def saveali(nres1, src, dst, totscore, bestscore_container):
    """
    Save alignment: copy src to dst, update bestscore.

    bestscore_container is a list [bestscore] to allow mutation.
    """
    for i in range(1, nres1 + 1):
        dst[i] = src[i]
    bestscore_container[0] = totscore


def initialize_lean(nres1, nres2):
    """
    Initialize all data structures for lean_mc.

    Returns dict of all state arrays.
    """
    # All arrays are 1-based. We allocate extra space.
    sz1 = nres1 + 2
    sz2 = nres2 + 2

    ali1 = np.zeros(sz1, dtype=np.int16)
    ali2 = np.zeros(sz2, dtype=np.int16)
    fragali1 = np.zeros(sz1, dtype=np.int16)
    fragali2 = np.zeros(sz2, dtype=np.int16)
    cnt = np.zeros(sz1, dtype=np.int32)

    prev1 = np.zeros(sz1, dtype=np.int32)
    next1 = np.full(sz1, nres1 + 1, dtype=np.int32)
    nextres = np.full(sz1, nres1 + 1, dtype=np.int32)
    next2 = np.full(sz2, nres2 + 1, dtype=np.int32)

    # ltest[i,j] for i=1..nres1, j=1..nres2
    ltest = np.zeros((sz1, sz2), dtype=np.bool_)

    # rescore and dscore matrices
    rescore = np.zeros((sz1, sz2), dtype=np.float32)
    dscore = np.zeros((sz1, sz2), dtype=np.float32)

    ncand = 0
    candij = np.zeros((3, MAXPAIR * 4 + 1), dtype=np.int32)  # candij[1..2, 1..ncand]

    return {
        'ali1': ali1, 'ali2': ali2,
        'fragali1': fragali1, 'fragali2': fragali2,
        'cnt': cnt,
        'prev1': prev1, 'next1': next1,
        'nextres': nextres, 'next2': next2,
        'ltest': ltest,
        'rescore': rescore, 'dscore': dscore,
        'ncand': ncand, 'candij': candij,
    }


def addscore(i, j, ali1, ali2, nextres, nres1, d1, d2):
    """
    Returns the marginal score of adding pair (i,j) relative to ali1.

    Translated from clean.f addscore().
    All arithmetic in float32 to match Fortran real precision.
    """
    zero = np.int16(0)
    dx = np.float32(scorefun(zero, zero))
    iold = ali2[j]
    jold = ali1[i]
    k = 0
    while nextres[k] <= nres1:
        k = nextres[k]
        x_val = ali1[k]
        if k != iold and x_val != jold:
            a = d1[i, k]
            b = d2[j, x_val]
            s = scorefun(a, b)
            # Fortran: dx=dx+s+s → (dx+s)+s (left-to-right)
            dx = np.float32(np.float32(dx + s) + s)
    return dx


def appendcand_lean(i, j, blocksize, state, nres1, d1, d2, lsave=True):
    """
    Add residue pair candidates to candidate list.

    Translated from lean.f appendcand_lean().
    """
    if not lsave:
        return

    ali1 = state['ali1']
    ali2 = state['ali2']
    ltest = state['ltest']
    dscore = state['dscore']
    candij = state['candij']
    nextres = state['nextres']

    for l in range(blocksize):
        j1 = abs(j + l)
        if not ltest[i + l, j1]:
            state['ncand'] += 1
            nc = state['ncand']
            candij[1, nc] = i + l
            candij[2, nc] = j1
            ltest[i + l, j1] = True

            if ali1[i + l] != j1:
                dscore[i + l, j1] = addscore(
                    i + l, j1, ali1, ali2, nextres, nres1, d1, d2)


def getrescore(i1, i2, d1, d2, ali1, nextres, nres1):
    """
    Calculate rescore(i1,i2) on the fly.

    Translated from clean.f getrescore().
    All arithmetic in float32 to match Fortran real precision.
    """
    zero = np.int16(0)
    dx = np.float32(scorefun(zero, zero))
    k = 0
    while nextres[k] <= nres1:
        k = nextres[k]
        x_val = ali1[k]
        if k != i1 and x_val != i2:
            s = scorefun(d1[i1, k], d2[i2, x_val])
            # Fortran: dx=dx+s+s → (dx+s)+s
            dx = np.float32(np.float32(dx + s) + s)
    return dx


def checktopo(i, j, blocksize, prev1, next1, fragali1, nres1):
    """
    Topology check: delete genes that violate sequential ordering.

    Translated from clean.f checktopo().

    Returns:
        genedel: list of gene positions to delete
    """
    genedel = []

    # Check forward: genes after (i+blocksize-1) with fragali < j-blocksize+1
    k = min(nres1 + 1, i + blocksize - 1)
    while next1[k] <= nres1:
        k = next1[k]
        l = fragali1[k]
        if l <= j - blocksize:
            genedel.append(k)
        else:
            break

    # Check backward: genes before (i-blocksize+1) with fragali >= j+blocksize
    k = max(0, i - blocksize + 1)
    while prev1[k] > 0:
        k = prev1[k]
        l = fragali1[k]
        if l >= j + blocksize:
            genedel.append(k)
        else:
            break

    return genedel


def geneoverlap(i, j, fragali1, fragali2, next1, next2,
                blocksize, nres1, nres2, existing_genedel):
    """
    Find genes that overlap with gene (i,j) in i- or j-direction.

    Translated from clean.f geneoverlap().

    Returns:
        additional genedel entries
    """
    genedel = list(existing_genedel)

    # Check i direction
    k = next1[max(0, i - blocksize)]
    while k <= min(nres1, i + blocksize - 1):
        l = fragali1[k]
        if k - l != i - j:
            genedel.append(k)
        k = next1[k]

    # Check j direction (parallel only — j > 0 in DALICON lean_mc)
    if j > 0:
        # parallel - parallel
        l = next2[max(0, j - blocksize)]
        while l <= min(nres2, j + blocksize - 1):
            k2 = fragali2[l]
            if k2 > 0 and (k2 - l != i - j) and (k2 >= i + blocksize or k2 <= i - blocksize):
                genedel.append(k2)
            l = next2[l]

        # parallel - antiparallel: always clash
        l = next2[max(0, j - 1)]
        while l <= min(nres2, j + blocksize + blocksize - 1):
            k2 = fragali2[l]
            if k2 < 0 and (-k2 >= i + blocksize or -k2 <= i - blocksize):
                genedel.append(-k2)
            l = next2[l]

    return genedel


def gowithgenes(i1, j1, genedel, ali1, cnt, blocksize, laddition):
    """
    Find residues whose cnt would drop to 0 after gene deletion.

    Translated from clean.f gowithgenes().

    Returns:
        rem: list of (res_i, res_j) tuples to remove
    """
    # Copy cnt to tmp for genes being deleted
    tmp = {}
    for gene_i in genedel:
        for l in range(blocksize):
            k = abs(gene_i + l)
            tmp[k] = int(cnt[k])

    # Subtract counts and find residues to remove
    rem = []
    for gene_i in genedel:
        for l in range(blocksize):
            k = abs(gene_i + l)
            tmp[k] -= 1
            if tmp[k] == 0:
                if laddition:
                    # Skip if this residue will be re-added by new gene (i1,j1)
                    m = k - i1
                    if 0 <= m <= blocksize - 1 and ali1[k] == abs(j1 + m):
                        continue
                rem.append((k, int(ali1[k])))

    return rem


def testaddition(i, j, blocksize, state, nres1, nres2, d1, d2, ltop=True, lsave=True):
    """
    Test adding gene (i,j): compute score change dx, residues to remove/add.

    Translated from clean.f testaddition().

    Returns:
        dx: score change
        nrem, rem: residues to delete [(i,j), ...]
        nnew, new: residues to add [(i,j), ...]
        genedel: genes to delete [i, ...]
    """
    prev1 = state['prev1']
    next1 = state['next1']
    ali1 = state['ali1']
    cnt = state['cnt']
    rescore = state['rescore']
    dscore = state['dscore']
    fragali1 = state['fragali1']
    fragali2 = state['fragali2']
    nextres = state['nextres']
    next2 = state['next2']

    genedel = []

    # Topology check
    if ltop:
        genedel = checktopo(i, j, blocksize, prev1, next1, fragali1, nres1)

    # Gene overlap
    genedel = geneoverlap(i, j, fragali1, fragali2, next1, next2,
                          blocksize, nres1, nres2, genedel)

    # Which residues go with deleted genes?
    rem = gowithgenes(i, j, genedel, ali1, cnt, blocksize, laddition=True)

    # New residues
    new = []
    for l in range(blocksize):
        if ali1[i + l] != abs(j + l):
            new.append((i + l, abs(j + l)))

    # Calculate dx — all float32 to match Fortran real precision
    dx = np.float32(0.0)

    # Add new
    for inew_idx, (i1, i2) in enumerate(new):
        if lsave:
            dx = np.float32(dx + dscore[i1, i2])
        else:
            dx = np.float32(dx + getrescore(i1, i2, d1, d2, ali1, nextres, nres1))

        # Add new x new
        for k in range(inew_idx):
            s = scorefun(d1[i1, new[k][0]], d2[i2, new[k][1]])
            # Fortran: dx=dx+s+s → (dx+s)+s
            dx = np.float32(np.float32(dx + s) + s)

        # Subtract new x rem
        for (remi, remj) in rem:
            if i1 != remi and i2 != remj:
                s = scorefun(d1[i1, remi], d2[i2, remj])
                # Fortran: dx=dx-s-s → (dx-s)-s
                dx = np.float32(np.float32(dx - s) - s)

    # Subtract rem x rem
    for k_idx, (remi, remj) in enumerate(rem):
        if lsave:
            dx = np.float32(dx - rescore[remi, remj])
        else:
            dx = np.float32(dx - getrescore(remi, remj, d1, d2, ali1, nextres, nres1))
        for l_idx in range(k_idx):
            s = scorefun(d1[remi, rem[l_idx][0]], d2[remj, rem[l_idx][1]])
            # Fortran: dx=dx+s+s → (dx+s)+s
            dx = np.float32(np.float32(dx + s) + s)

    return dx, rem, new, genedel


def testdeletion(genedel, blocksize, state, nres1, d1, d2, lsave=True):
    """
    Test deleting genes: compute score change dx.

    Translated from clean.f testdeletion().

    Returns:
        dx: score change (positive means improvement from deletion)
        rem: residues to delete [(i,j), ...]
    """
    ali1 = state['ali1']
    cnt = state['cnt']
    rescore = state['rescore']
    nextres = state['nextres']

    all_rem = []
    dx = np.float32(0.0)

    for gene_i in genedel:
        j = ali1[gene_i]
        rem = gowithgenes(gene_i, j, [gene_i], ali1, cnt, blocksize, laddition=False)

        for l_idx, (i1, i2) in enumerate(rem):
            if lsave:
                dx = np.float32(dx - rescore[i1, i2])
            else:
                dx = np.float32(dx - getrescore(i1, i2, d1, d2, ali1, nextres, nres1))

            for k_idx in range(l_idx):
                j1, j2 = rem[k_idx]
                s = scorefun(d1[i1, j1], d2[i2, j2])
                # Fortran: dx=dx+s+s → (dx+s)+s
                dx = np.float32(np.float32(dx + s) + s)

        all_rem.extend(rem)

    return dx, all_rem


def dodeletions_lean(genedel, rem, state, blocksize, d1, d2, lsave=True):
    """
    Execute gene and residue deletions.

    Translated from lean.f dodeletions_lean().
    """
    fragali1 = state['fragali1']
    fragali2 = state['fragali2']
    ali1 = state['ali1']
    ali2 = state['ali2']
    cnt = state['cnt']
    prev1 = state['prev1']
    next1 = state['next1']
    next2 = state['next2']
    nextres = state['nextres']
    rescore = state['rescore']
    dscore = state['dscore']
    candij = state['candij']
    ncand = state['ncand']

    # Delete genes
    for gene_i in genedel:
        j = int(fragali1[gene_i])
        fragali1[gene_i] = 0
        fragali2[abs(j)] = 0
        for l in range(blocksize):
            cnt[gene_i + l] -= 1

        # Update prev1, next1
        ip = prev1[gene_i]
        in_ = next1[gene_i]
        for ix in range(gene_i, in_ + 1):
            if ix < len(prev1):
                prev1[ix] = ip
        for ix in range(ip, gene_i + 1):
            if ix >= 0:
                next1[ix] = in_

        # Update next2
        j1 = abs(j)
        in_ = next2[j1]
        ix = j1
        if ix > 0:
            while ix > 0 and next2[ix - 1] == j1:
                ix -= 1
                next2[ix] = in_
                if ix <= 0:
                    break

    # Delete residues
    ndel = 0
    for (ri, rj) in rem:
        # Update nextres
        in_ = nextres[ri]
        ix = ri
        while ix > 0 and nextres[ix - 1] == ri:
            ix -= 1
            nextres[ix] = in_
            if ix <= 0:
                break

        if lsave:
            # Update dscore/rescore — float32 arithmetic
            for k in range(1, ncand + 1):
                cix = candij[1, k]
                cjx = candij[2, k]
                ddx = np.float32(0.0)
                if ri != cix and rj != cjx:
                    s = np.float32(-scorefun(d1[cix, ri], d2[cjx, rj]))
                    ddx = np.float32(s + s)
                if ali1[cix] == cjx:
                    rescore[cix, cjx] = np.float32(rescore[cix, cjx] + ddx)
                else:
                    dscore[cix, cjx] = np.float32(dscore[cix, cjx] + ddx)
            dscore[ri, rj] = rescore[ri, rj]
            rescore[ri, rj] = np.float32(0.0)

        ali1[ri] = 0
        ali2[rj] = 0
        ndel += 1

    return ndel


def doadditions_lean(i, j, new, state, blocksize, d1, d2, nres1, lsave=True):
    """
    Execute gene and residue additions.

    Translated from lean.f doadditions_lean().
    """
    fragali1 = state['fragali1']
    fragali2 = state['fragali2']
    ali1 = state['ali1']
    ali2 = state['ali2']
    cnt = state['cnt']
    prev1 = state['prev1']
    next1 = state['next1']
    next2 = state['next2']
    nextres = state['nextres']
    rescore = state['rescore']
    dscore = state['dscore']
    candij = state['candij']
    ncand = state['ncand']

    # Add gene
    fragali1[i] = j
    if j < 0:
        fragali2[-j] = -i
    elif j > 0:
        fragali2[j] = i
    for l in range(blocksize):
        cnt[i + l] += 1

    # Update prev1, next1
    ip = prev1[i]
    in_ = next1[i]
    for ix in range(ip, i):
        if ix >= 0:
            next1[ix] = i
    for ix in range(i + 1, in_ + 1):
        if ix < len(prev1):
            prev1[ix] = i

    # Update next2
    j1 = abs(j)
    in_ = next2[j1]
    ix = j1
    if ix > 0:
        while ix > 0 and next2[ix - 1] == in_:
            ix -= 1
            next2[ix] = j1
            if ix <= 0:
                break

    # Add residues
    nacc = 0
    for (i1, j1_val) in new:
        # Update nextres
        in_ = nextres[i1]
        ix = i1
        while ix > 0 and nextres[ix - 1] == in_:
            ix -= 1
            nextres[ix] = i1
            if ix <= 0:
                break

        if lsave:
            # Update dscore/rescore — float32 arithmetic
            for k in range(1, ncand + 1):
                cix = candij[1, k]
                cjx = candij[2, k]
                ddx = np.float32(0.0)
                if i1 != cix and j1_val != cjx:
                    s = scorefun(d1[cix, i1], d2[cjx, j1_val])
                    ddx = np.float32(s + s)
                if ali1[cix] == cjx:
                    rescore[cix, cjx] = np.float32(rescore[cix, cjx] + ddx)
                else:
                    dscore[cix, cjx] = np.float32(dscore[cix, cjx] + ddx)
            rescore[i1, j1_val] = dscore[i1, j1_val]
            dscore[i1, j1_val] = np.float32(0.0)

        # Handle potential conflicts
        iold = ali2[j1_val]
        jold = ali1[i1]
        if iold != 0:
            ali1[iold] = 0
        if jold != 0:
            ali2[jold] = 0
        ali1[i1] = j1_val
        ali2[j1_val] = i1
        nacc += 1

    return nacc


def changeconfig_lean(genedel, rem, i, j, new, state, blocksize,
                      d1, d2, nres1, dx, totscore, bestscore,
                      bestali1, outfragali, lsave=True):
    """
    Apply changes: deletions then additions, update totscore, save best.

    Translated from lean.f changeconfig_lean().

    Returns:
        ndel, nacc: counts
        totscore: updated total score
        bestscore: updated best score
    """
    ndel = dodeletions_lean(genedel, rem, state, blocksize, d1, d2, lsave)
    nacc = 0
    if i != 0:
        nacc = doadditions_lean(i, j, new, state, blocksize, d1, d2, nres1, lsave)

    totscore = np.float32(totscore + dx)

    if totscore > bestscore[0]:
        saveali(nres1, state['fragali1'], outfragali, totscore, bestscore)
        saveali(nres1, state['ali1'], bestali1, totscore, bestscore)

    return ndel, nacc, totscore


def lean_mc(imax0, blocksize, nres1, nres2, d1, d2, lsave,
            preali1, itrimx, ntetra, tetrapool, rng):
    """
    Main lean Monte Carlo optimization.

    Translated from lean.f lean_mc().

    Args:
        imax0: max Monte Carlo steps
        blocksize: alignment block size (4 for DALICON)
        nres1, nres2: sequence lengths
        d1, d2: distance matrices (1-based)
        lsave: use rescore/dscore caching
        preali1: prealignment (1-based, int16)
        itrimx: trim frequency
        ntetra: number of tetrapeptide candidates
        tetrapool: list of (i,j) tuples
        rng: RandomState instance

    Returns:
        score: best score / 100
        bestali1: best alignment array (1-based)
    """
    imax = imax0

    # Allocate result arrays
    sz1 = nres1 + 2
    bestali1 = np.zeros(sz1, dtype=np.int16)
    outfragali = np.zeros(sz1, dtype=np.int16)
    bestscore = [np.float32(0.0)]  # mutable container, float32

    # Initialize
    state = initialize_lean(nres1, nres2)

    # Load candidates
    for ix in range(ntetra):
        i, j = tetrapool[ix]
        appendcand_lean(i, j, blocksize, state, nres1, d1, d2, lsave)

    # Convert preali1 to initfragali
    initfragali = np.zeros(sz1, dtype=np.int16)
    for i in range(1, min(nres1 + 1, len(preali1))):
        initfragali[i] = 0

    for i in range(1, min(nres1 - 2, len(preali1))):
        j = int(preali1[i])
        if j > 0 and j <= nres2 - 3:
            # At least one candidate per segment, more on long ones
            if i > 1:
                if preali1[i - 1] == 0:
                    initfragali[i] = j
            if i + 3 < len(preali1) and preali1[i + 3] == j + 3:
                initfragali[i] = j

    # Load initial alignment
    lenali = 0
    totscore = np.float32(0.0)
    for i in range(1, nres1 + 1):
        j = int(initfragali[i])
        if j != 0:
            dx, rem, new, genedel = testaddition(
                i, j, blocksize, state, nres1, nres2, d1, d2,
                ltop=True, lsave=lsave)
            # Force accept (aconitase fix)
            bestscore[0] = np.float32(0.0)
            ndel, nacc, totscore = changeconfig_lean(
                genedel, rem, i, j, new, state, blocksize,
                d1, d2, nres1, dx, totscore, bestscore,
                bestali1, outfragali, lsave)

    # Monte Carlo steps — all float32 to match Fortran real precision
    imp = 0
    istep = 0
    impscore = np.float32(0.0)
    alfa = np.float32(50.0)
    delta = np.float32(0.0)

    # Allocate sort arrays
    di = np.zeros(ntetra + 1, dtype=np.int32)
    ve = np.zeros(ntetra + 1, dtype=np.int32)

    while istep < imax:
        istep += 1

        # Heat shock
        if istep - imp == 1:
            alfa = np.float32(50.0)
            delta = np.float32(0.0)
        elif istep - imp == 2:
            alfa = np.float32(np.float32(1.0) / np.sqrt(max(np.float32(10000.0), totscore)))
            delta = np.float32(alfa / np.float32(1000.0))

        ndel_total = 0
        nacc_total = 0
        alfa = np.float32(alfa + delta)

        # Adding phase (WARM) — test in randomized order
        for ix in range(1, ntetra + 1):
            di[ix] = ix
            ve[ix] = _nint(float(np.float32(100.0) * rng.ran()))

        _index_qsort(di, ve, 1, ntetra)

        for q in range(1, ntetra + 1):
            ix = di[q]
            i, j = tetrapool[ix - 1]  # tetrapool is 0-based list
            if state['fragali1'][i] != j:
                dx, rem, new, genedel = testaddition(
                    i, j, blocksize, state, nres1, nres2, d1, d2,
                    ltop=True, lsave=lsave)
                p = getp(np.float32(alfa * dx))
                if rng.ran() < p:
                    ndel, nacc, totscore = changeconfig_lean(
                        genedel, rem, i, j, new, state, blocksize,
                        d1, d2, nres1, dx, totscore, bestscore,
                        bestali1, outfragali, lsave)
                    ndel_total += ndel
                    nacc_total += nacc
                    if float(totscore) > float(impscore):
                        imp = istep
                        impscore = np.float32(totscore)

        # Trimming phase (COLD)
        if istep > 0:
            i = 0
            while state['next1'][i] <= nres1:
                i = state['next1'][i]
                j = int(state['fragali1'][i])
                # Trim from ends only
                if cnt_at_end(state['cnt'], i, blocksize):
                    genedel_trim = [i]
                    dx_trim, rem_trim = testdeletion(
                        genedel_trim, blocksize, state, nres1, d1, d2, lsave)
                    # Only accept improvements
                    if float(dx_trim) > 0.0:
                        ndel, nacc, totscore = changeconfig_lean(
                            genedel_trim, rem_trim, 0, 0, [], state, blocksize,
                            d1, d2, nres1, dx_trim, totscore, bestscore,
                            bestali1, outfragali, lsave)
                        ndel_total += ndel
                        nacc_total += nacc
                        if float(totscore) > float(np.float32(impscore + np.float32(1.0))):
                            imp = istep
                            impscore = np.float32(totscore)

        # Early return on hopeless cases
        ncand = state['ncand']
        if nacc_total > ncand:
            return float(bestscore[0]) / 100.0, bestali1
        if ndel_total > ncand:
            return float(bestscore[0]) / 100.0, bestali1

        # Quit if no improvement for 20 steps
        if istep - imp >= 20:
            break

        # Accelerate cooling if many acceptances
        if ntetra > 0 and nacc_total / float(ntetra) > 0.1:
            alfa = np.float32(np.float32(10.0) * alfa)

    return float(bestscore[0]) / 100.0, bestali1


def cnt_at_end(cnt, i, blocksize):
    """Check if gene at position i has cnt==1 at either end."""
    return cnt[i] == 1 or cnt[i + blocksize - 1] == 1
