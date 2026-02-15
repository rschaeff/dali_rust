# DSSP Entry Points in DaliLite.v5

## Overview

DSSP (Dictionary of Secondary Structure of Proteins) serves as a preprocessor
that reduces a full-atom PDB structure to an SSE-annotated CA trace. This
document traces how DSSP data enters and flows through the Dali comparison
pipeline.

## Pipeline: PDB to .dat

```
PDB file
  |
  v
dsspcmbi (src/dsspcmbi/)          -- bundled C implementation of Kabsch & Sander (1983)
  |                                   125K lines of auto-translated Pascal->C
  v
.dssp file                         -- standard DSSP output: per-residue 8-state SS, CA coords,
  |                                   H-bond partners, accessibility, backbone angles
  |
  +---> puu (src/puu.f)            -- reads DSSP via readdssp1() -> getdssp() (src/util.f:155)
  |       extracts: struc, ca, seq, H-bond partners (bp0)
  |       builds domain decomposition tree -> units.puu
  |
  +---> puutos (src/puutos.f)      -- reads DSSP again via readdssp1()
          calls getsecstr() to reduce 8-state -> 3-state
          writes .dat file with SSE annotation
```

Perl orchestration in `mpidali.pm:prepare_brkfile()` (line 591):
```perl
$cmd = "$DSSP_EXE $query_pdbfile $query_dsspfile";    # line 624
$cmd = "$PUU_EXE $shortcode[$i] $query_pdbfile $query_dsspfile";  # line 699
$cmd = "$PUUTOS_EXE units.puu $query_dsspfile $DALIDATDIR/";      # line 702
```

## The Critical Transformation: getsecstr()

`puutos.f:577` — converts DSSP 8-state to simplified H/E/L and filters by length.

### 8-state to 3-state mapping

```fortran
character*5 fromstring,tostring
data fromstring/'GIHE*'/
data tostring  /'HHHEL'/
```

| DSSP code | Meaning              | Mapped to |
|-----------|----------------------|-----------|
| G         | 3-10 helix           | H (helix) |
| I         | pi helix             | H (helix) |
| H         | alpha helix          | H (helix) |
| E         | extended strand      | E (strand)|
| * (else)  | everything else      | L (loop)  |

B (bridge), T (turn), S (bend), blank (coil) all map to L.

### Minimum length filtering

- Loop segments: **deleted entirely** (only H and E survive)
- Helices shorter than **6 residues**: excluded
- Strands shorter than **8 residues**: excluded (note: the call is `getsecstr(nres,struc,nseg,secstr,segmentrange,6,8)` where args are lene=6, lenh=8 — naming is counterintuitive but lenh is the helix threshold)
- Borderline SSEs (above absolute minimum but below threshold): extended into adjacent loop regions via `growleft()`/`growrite()`

### Output

- `secstr(nseg)` — array of 'H' or 'E', one per SSE
- `segmentrange(2,nseg)` — start/end residue index per SSE
- `na` — count of helices
- `nb` — count of strands

## .dat File Format

```
>>>> 101mA  154    5    5    0  HHHHH
         1         5        36         5        36         0     <- segmentrange + checkrange + checkx
         2        60        77        60        77         0
         ...
     CA coordinates (10f8.1)
>>>> 101mA  154    7                                             <- domain tree
     domain nodes (idom, node_type, children, segment lists)
>>>> 101mA  154    7                                             <- domain tree (compact form)
-dssp     "LLLHHHHHHHH..."                                      <- full per-residue DSSP string
-sequence "VLSEGEWQLVL..."                                      <- amino acid sequence
-compnd   "MOLECULE: MYOGLOBIN"                                 <- PDB metadata
-acc  ...                                                       <- per-residue accessibility + CA coords
```

Header format: `>>>> code chainid nres nseg na nb secstr(1..nseg)`

## How SSE Data Drives Comparison

### WOLF/SOAP routing (mpidali.pm)

```perl
my $MINSSE = 3;
($nres, $nsse) = />>>>\s+\S+\s+(\d+)\s+(\d+)/;   # parse .dat header
if ($nsse < $MINSSE) { &wolf("SOAP", ...); }       # Needleman-Wunsch on CA distances
else                 { &wolf("WOLF", ...); }        # SSE vector hashing
```

Examples: 1bbaA (nseg=1, 1 helix) and 1pptA (nseg=1) route to SOAP.
101mA (nseg=5, all helices) and 1a87A (nseg=17) route to WOLF.

### WOLF comparison (wolf_original.f:compare, line 31)

SSE type is a **hard filter** — only same-type pairs are considered:

```fortran
if(protein_secstr(cur_a).ne.atype) goto 19   ! skip
if(protein_secstr(cur_c).ne.ctype) goto 19   ! only H<->H, E<->E
```

SSE midpoints and direction vectors (computed from CA coordinates of segment
residues) are hashed into a 3D spatial grid for rapid neighbor lookup.

### PARSI viability check (parsi-admin.f:getstructure1, line 188)

```fortran
if ((nseg.le.2).or.(ndom.lt.nseg).or.(na1*na2+nb1*nb2.eq.0)) then
    flag=1   ! skip — insufficient SSEs or no matching types
    return
end if
```

Requires: nseg > 2, valid domain tree, AND at least one shared SSE type
(both proteins have helices, or both have strands).

### FILTER95 (comparemodules.f)

Uses `secstr` and `checkrange` for domain-aware Z-score computation with
SSE-sensitive boundary detection.

## What Dali Actually Consumes from DSSP

| Field | Extracted by | Used by |
|-------|-------------|---------|
| `secstr(nseg)` — H/E type per SSE | `getsecstr()` | WOLF type filter, PARSI viability, FILTER95 boundaries |
| `segmentrange(2,nseg)` — residue ranges | `getsecstr()` | WOLF vectors, PARSI search space definition |
| `ca(3,nres)` — CA coordinates | `getdssp()` column 83+ | All stages (distance matrices, superposition, scoring) |
| `na`, `nb` — helix/strand counts | `countchr()` on secstr | WOLF/SOAP routing, PARSI viability |
| H-bond partners (bp0) | `getdssp()` | `puu` domain decomposition only, NOT comparison |
| Accessibility, angles | `getdssp()` | Stored in .dat but unused by comparison algorithms |

## Reimplementation Notes

### Not a priority for initial reimplementation

The import pipeline (PDB -> DSSP -> .dat) is upstream of the comparison
algorithm. For validation purposes, we can use the existing .dat files
directly — the comparison algorithm reads .dat, not DSSP.

### Drop-in replacement is straightforward

The comparison algorithm only needs the 3-state reduction produced by
`getsecstr()`. Any DSSP implementation that produces the standard 8-state
assignment would work:

- **mkdssp** (modern C++ DSSP, maintained by PDB-REDO group)
- **BioPython** `Bio.PDB.DSSP` (wrapper around mkdssp)
- **pydssp** (pure Python reimplementation)
- **biotite** `annotate_sse()` (pure Python, different algorithm)

The `getsecstr()` logic itself is ~100 lines of straightforward Fortran
that maps 8-state to 3-state and filters by length.

### Future: CIF format support

PDB format structure representations are becoming obsolete. The PDB archive
has fully transitioned to mmCIF as the canonical format, and many newer
structures (large complexes, structures from AlphaFold/ESMFold) are only
available in CIF.

A future reimplementation of the import pathway should:

1. **Read mmCIF** (PDBx/mmCIF) as the primary input format, with legacy PDB
   as a fallback. Libraries: `gemmi` (C++/Python, fast), `BioPython`
   `MMCIFParser`, `atomium`, or direct parsing of the `_atom_site` category.

2. **Compute secondary structure from coordinates** rather than requiring a
   separate DSSP binary. Options:
   - `mkdssp` can read mmCIF directly (v4+)
   - `pydssp` operates on coordinate arrays
   - `biotite` provides `annotate_sse()` from coordinate arrays
   - Reimplement Kabsch-Sander H-bond criteria directly (the algorithm is
     well-defined: E_hb = 0.084 * {1/rON + 1/rCH - 1/rOH - 1/rCN} * 332 kcal/mol,
     H-bond if E < -0.5 kcal/mol)

3. **Preserve the `getsecstr()` reduction** exactly — the 8-to-3 state
   mapping and minimum length thresholds (6 residues for helices, 8 for
   strands) are algorithmic parameters that affect downstream results.

4. **Decouple the .dat format** from the import step. The .dat format itself
   is a reasonable internal representation (SSE list + CA trace + domain tree),
   but could be replaced with a more modern serialization (e.g., HDF5,
   MessagePack, or even a Python dataclass serialized with pickle/JSON).

### Key source files

| File | Role |
|------|------|
| `src/dsspcmbi/DsspCMBI.c` | Bundled DSSP implementation (125K lines C) |
| `src/util.f:getdssp()` (line 155) | DSSP file parser |
| `src/puutos.f:getsecstr()` (line 577) | 8-state to 3-state reduction + length filtering |
| `src/puutos.f:writeit()` (line 76) | Writes .dat file |
| `src/puu.f:puu1()` | Domain decomposition from DSSP |
| `src/wolf_original.f:readproteindata()` (line 374) | Reads .dat SSE fields |
| `src/wolf_original.f:compare()` (line 31) | Uses secstr as hard match filter |
| `src/parsi-admin.f:getstructure1()` (line 150) | SSE viability check for PARSI |
| `bin/mpidali.pm:prepare_brkfile()` (line 591) | Perl orchestration of import |
| `bin/mpidali.pm:get_size_nsse()` (line 822) | Reads nseg from .dat for WOLF/SOAP routing |
