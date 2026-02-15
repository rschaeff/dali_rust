"""
3D spatial hash grid for SSE descriptor lookup.

Replaces Fortran linked-list grid with Python defaultdict.
Grid cell size = 2.0 Å, bounds = -20..+20 in each dimension.
"""

from collections import defaultdict
import numpy as np


def fung(x):
    """Quantize coordinate to grid index. Grid spacing = 2 Å."""
    # Fortran nint: round half to nearest even → use np.rint for compatibility
    return int(np.rint(x / 2.0))


def lgrid(gx, gy, gz, maxgrid=20):
    """Check if grid indices are within bounds."""
    return (abs(gx) <= maxgrid and abs(gy) <= maxgrid and abs(gz) <= maxgrid)


class SpatialHashGrid:
    """
    3D spatial hash grid storing SSE descriptors.

    Each entry stores:
      - aseg: first SSE of the pair that defined the canonical frame (0-based)
      - bseg: second SSE of the pair (0-based)
      - cseg: the "other" SSE being hashed (0-based)
      - link_from: (3,) transformed N-terminal endpoint
      - link_to: (3,) transformed C-terminal endpoint
    """

    def __init__(self):
        self.grid = defaultdict(list)

    def clear(self):
        self.grid.clear()

    def boxit(self, nseg, x, aseg, bseg):
        """
        Insert all SSE descriptors (except aseg) into the grid.

        Args:
            nseg: number of SSEs
            x: (3, 3 + 2*nseg) transformed coordinate array (after twist)
            aseg: SSE index to exclude (0-based)
            bseg: paired SSE index (0-based)
        """
        for iseg in range(nseg):
            if iseg == aseg:
                continue

            # Transformed endpoints
            link_from = x[:, 3 + iseg].copy()
            link_to = x[:, 3 + nseg + iseg].copy()

            # Grid position = midpoint of transformed endpoints
            mid = (link_from + link_to) / 2.0
            gx = fung(mid[0])
            gy = fung(mid[1])
            gz = fung(mid[2])

            if lgrid(gx, gy, gz):
                entry = (aseg, bseg, iseg, link_from, link_to)
                self.grid[(gx, gy, gz)].append(entry)

    def lookup(self, gx, gy, gz, radius=2):
        """
        Return all entries in grid cells within ±radius of (gx, gy, gz).

        Yields: (aseg, bseg, cseg, link_from, link_to) tuples
        """
        for dx in range(-radius, radius + 1):
            for dy in range(-radius, radius + 1):
                for dz in range(-radius, radius + 1):
                    key = (gx + dx, gy + dy, gz + dz)
                    if key in self.grid:
                        yield from self.grid[key]
