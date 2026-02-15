"""Constants from parsizes.for."""

MAXRES1 = 4000      # max residues in protein 1 (compressed)
MAXRES2 = 4000      # max residues in protein 2
MAXRES0 = 400       # maxres2/10 — max candidates per segment
EXDIM = 1000000     # size of score table ex()
MAXSTACK = 10000    # stack size (in integers)
MAXSTACK2 = 10000   # link array size
MAXSEG = 200        # max SSE segments
MAXDOM = 400        # max domain tree nodes (maxseg*2)
INFINIT = 10000000  # sentinel for "infinity"
NUL = -99           # sentinel for "unaligned"
BOXDIM = 501        # link_pointer grid dimension
BL = 10             # block stepping size
STARTSIZE = 3       # minimum segments to align a domain node
LENE = 6            # min segment length for E (strand)
LENH = 8            # min segment length for H (helix)
