use ndarray::Array2;

/// Distance matrix scaled by 10 (used by DP module).
/// Computed as nint(dist * 10), stored as i16.
#[derive(Debug, Clone)]
pub struct DistMatScale10 {
    pub data: Array2<i16>,
    pub nres: usize,
}

/// Distance matrix scaled by 100 (used by FILTER95 module).
/// Computed as nint(dist * 100), stored as i16.
#[derive(Debug, Clone)]
pub struct DistMatScale100 {
    pub data: Array2<i16>,
    pub nres: usize,
}

/// Distance matrix for DALICON (1-based indexing convention).
#[derive(Debug, Clone)]
pub struct DaliconDistMat {
    pub data: Array2<i16>,
    pub nres: usize,
}

/// Distance matrix for PARSI (capped at 1000).
#[derive(Debug, Clone)]
pub struct ParsiDistMat {
    pub data: Array2<i16>,
    pub nres: usize,
}
