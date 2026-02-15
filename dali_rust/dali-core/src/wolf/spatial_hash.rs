use std::collections::HashMap;
use ndarray::Array2;

/// Grid spacing in Angstroms.
const GRID_SPACING: f64 = 2.0;
/// Maximum grid index in each dimension.
const MAX_GRID: i32 = 20;

/// Quantize coordinate to grid index. Grid spacing = 2 Å.
#[inline]
pub fn fung(x: f64) -> i32 {
    // Fortran nint / np.rint: round half to even
    // For grid quantization, this matches the Python implementation
    (x / GRID_SPACING).round() as i32
}

/// Check if grid indices are within bounds.
#[inline]
pub fn lgrid(gx: i32, gy: i32, gz: i32) -> bool {
    gx.abs() <= MAX_GRID && gy.abs() <= MAX_GRID && gz.abs() <= MAX_GRID
}

/// Entry stored in the spatial hash grid.
#[derive(Clone)]
pub struct GridEntry {
    pub aseg: usize,       // first SSE of the pair that defined the canonical frame
    pub bseg: usize,       // second SSE of the pair
    pub cseg: usize,       // the "other" SSE being hashed
    pub link_from: [f64; 3], // transformed N-terminal endpoint
    pub link_to: [f64; 3],   // transformed C-terminal endpoint
}

/// 3D spatial hash grid storing SSE descriptors.
pub struct SpatialHashGrid {
    grid: HashMap<(i32, i32, i32), Vec<GridEntry>>,
}

impl SpatialHashGrid {
    pub fn new() -> Self {
        SpatialHashGrid {
            grid: HashMap::new(),
        }
    }

    pub fn clear(&mut self) {
        self.grid.clear();
    }

    /// Insert all SSE descriptors (except aseg) into the grid.
    ///
    /// # Arguments
    /// * `nseg` - number of SSEs
    /// * `x` - (3, 3 + 2*nseg) transformed coordinate array (after twist)
    /// * `aseg` - SSE index to exclude (0-based)
    /// * `bseg` - paired SSE index (0-based)
    pub fn boxit(&mut self, nseg: usize, x: &Array2<f64>, aseg: usize, bseg: usize) {
        for iseg in 0..nseg {
            if iseg == aseg {
                continue;
            }

            // Transformed endpoints
            let link_from = [
                x[[0, 3 + iseg]],
                x[[1, 3 + iseg]],
                x[[2, 3 + iseg]],
            ];
            let link_to = [
                x[[0, 3 + nseg + iseg]],
                x[[1, 3 + nseg + iseg]],
                x[[2, 3 + nseg + iseg]],
            ];

            // Grid position = midpoint of transformed endpoints
            let mid_x = (link_from[0] + link_to[0]) / 2.0;
            let mid_y = (link_from[1] + link_to[1]) / 2.0;
            let mid_z = (link_from[2] + link_to[2]) / 2.0;

            let gx = fung(mid_x);
            let gy = fung(mid_y);
            let gz = fung(mid_z);

            if lgrid(gx, gy, gz) {
                let entry = GridEntry {
                    aseg,
                    bseg,
                    cseg: iseg,
                    link_from,
                    link_to,
                };
                self.grid.entry((gx, gy, gz)).or_default().push(entry);
            }
        }
    }

    /// Iterate over all entries in grid cells within ±radius of (gx, gy, gz).
    pub fn lookup(&self, gx: i32, gy: i32, gz: i32, radius: i32) -> GridLookupIter<'_> {
        GridLookupIter {
            grid: &self.grid,
            gx,
            gy,
            gz,
            radius,
            dx: -radius,
            dy: -radius,
            dz: -radius,
            entry_idx: 0,
            current_entries: None,
        }
    }
}

impl Default for SpatialHashGrid {
    fn default() -> Self {
        Self::new()
    }
}

/// Iterator over grid lookup results.
pub struct GridLookupIter<'a> {
    grid: &'a HashMap<(i32, i32, i32), Vec<GridEntry>>,
    gx: i32,
    gy: i32,
    gz: i32,
    radius: i32,
    dx: i32,
    dy: i32,
    dz: i32,
    entry_idx: usize,
    current_entries: Option<&'a Vec<GridEntry>>,
}

impl<'a> Iterator for GridLookupIter<'a> {
    type Item = &'a GridEntry;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // Try to return next entry from current cell
            if let Some(entries) = self.current_entries {
                if self.entry_idx < entries.len() {
                    let entry = &entries[self.entry_idx];
                    self.entry_idx += 1;
                    return Some(entry);
                }
                self.current_entries = None;
            }

            // Advance to next cell
            if self.dx > self.radius {
                return None;
            }

            let key = (self.gx + self.dx, self.gy + self.dy, self.gz + self.dz);

            // Advance dz, dy, dx
            self.dz += 1;
            if self.dz > self.radius {
                self.dz = -self.radius;
                self.dy += 1;
                if self.dy > self.radius {
                    self.dy = -self.radius;
                    self.dx += 1;
                }
            }

            // Check if this cell has entries
            if let Some(entries) = self.grid.get(&key) {
                self.current_entries = Some(entries);
                self.entry_idx = 0;
            }
        }
    }
}
