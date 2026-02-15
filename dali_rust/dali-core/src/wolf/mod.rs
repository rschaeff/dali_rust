pub mod geometry;
pub mod spatial_hash;
pub mod compare;

use std::path::Path;
use ndarray::Array2;

use crate::io::dat;
use crate::types::Protein;
use crate::types::{AlignmentBlock, compress_blocks};
use crate::numerics::fitz::fitz;
use geometry::{vec_sse, twist, preparex, compute_neidist};
use spatial_hash::SpatialHashGrid;
use compare::compare;

/// Constants matching Fortran serialcompare.f
const RCUT_FITZ: f64 = 4.0;
const MAXITER: usize = 20;
const NEIBORCUTOFF: f64 = 12.0;

/// Pre-computed data for a protein's SSE geometry.
pub struct WolfProteinData {
    pub protein: Protein,
    pub midpoint: Array2<f64>,  // (3, nseg)
    pub direction: Array2<f64>, // (3, nseg)
    pub neidist: Array2<f64>,   // (nseg, nseg)
}

/// Result of a WOLF comparison.
#[derive(Debug, Clone)]
pub struct WolfResult {
    pub cd1: String,
    pub cd2: String,
    pub blocks: Vec<AlignmentBlock>,
}

impl WolfResult {
    /// Format as WOLFITZ output line matching Fortran free-format write.
    pub fn to_wolfitz_line(&self) -> String {
        let mut parts = Vec::new();
        parts.push(format!("WOLFITZ {:<5}{:<5}", self.cd1, self.cd2));
        let nblock = self.blocks.len();
        parts.push(format!("{:>12}", nblock));
        for b in &self.blocks {
            parts.push(format!("{:>12}{:>12}", b.l1, b.r1));
        }
        for b in &self.blocks {
            parts.push(format!("{:>12}{:>12}", b.l2, b.r2));
        }
        parts.join("")
    }
}

/// Read .dat file and compute SSE vectors and neighbor distances.
pub fn setup_protein<P: AsRef<Path>>(filepath: P) -> Result<WolfProteinData, dat::DatError> {
    let protein = dat::read_dat(filepath)?;
    let seg_starts: Vec<u32> = protein.segments.iter().map(|s| s.start).collect();
    let seg_ends: Vec<u32> = protein.segments.iter().map(|s| s.end).collect();

    let (midpoint, direction) = vec_sse(
        &protein.ca,
        protein.nres,
        protein.nseg,
        &seg_starts,
        &seg_ends,
    );
    let neidist = compute_neidist(protein.nseg, &midpoint);

    Ok(WolfProteinData {
        protein,
        midpoint,
        direction,
        neidist,
    })
}

/// Build spatial hash grid for a protein's SSE descriptors.
pub fn load_protein(
    grid: &mut SpatialHashGrid,
    data: &WolfProteinData,
) {
    grid.clear();
    let nseg = data.protein.nseg;
    for iseg in 0..nseg {
        for jseg in 0..nseg {
            if jseg != iseg && data.neidist[[iseg, jseg]] < NEIBORCUTOFF {
                let mut x = preparex(iseg, jseg, nseg, &data.midpoint, &data.direction);
                twist(&mut x, 3 + 2 * nseg);
                grid.boxit(nseg, &x, iseg, jseg);
            }
        }
    }
}

/// Compare two proteins using the WOLF SSE hashing algorithm.
///
/// # Arguments
/// * `cd1_data` - pre-computed data for protein 1
/// * `cd2_data` - pre-computed data for protein 2
/// * `grid` - pre-built SpatialHashGrid for cd1
pub fn wolf_compare(
    cd1_data: &WolfProteinData,
    cd2_data: &WolfProteinData,
    grid: &SpatialHashGrid,
) -> Option<WolfResult> {
    let prot1 = &cd1_data.protein;
    let prot2 = &cd2_data.protein;

    // Check both proteins have enough SSEs
    if prot1.nseg <= 2 || prot2.nseg <= 2 {
        return None;
    }

    // Compare: cd2 SSE pairs searched against cd1's grid
    // Note: Fortran passes nseg (cd2's) as protein_nseg
    let bp = compare(
        prot2.nseg,
        &cd2_data.midpoint,
        &cd2_data.direction,
        &prot2.secstr,
        prot2.nseg, // protein_nseg = cd2's nseg (Fortran behavior)
        &prot1.secstr,
        grid,
        NEIBORCUTOFF,
        &cd2_data.neidist,
    )?;

    if bp.iseg == bp.jseg || bp.useg == bp.vseg {
        return None;
    }

    // Build initial superposition for cd1 (using bestpair SSEs from cd1)
    let mut x_full = Array2::zeros((3, 3 + prot1.nres));
    {
        let x_ref = preparex(bp.useg, bp.vseg, 0, &cd1_data.midpoint, &cd1_data.direction);
        for k in 0..3 {
            for j in 0..3 {
                x_full[[k, j]] = x_ref[[k, j]];
            }
        }
    }
    for j in 0..prot1.nres {
        for k in 0..3 {
            x_full[[k, j + 3]] = prot1.ca[[k, j]];
        }
    }
    twist(&mut x_full, 3 + prot1.nres);

    // Extract transformed CA coords (skip 3 reference points)
    let mut xca = Array2::zeros((3, prot1.nres));
    for j in 0..prot1.nres {
        for k in 0..3 {
            xca[[k, j]] = x_full[[k, j + 3]];
        }
    }

    // Build initial superposition for cd2 (using bestpair SSEs from cd2)
    let mut y_full = Array2::zeros((3, 3 + prot2.nres));
    {
        let y_ref = preparex(bp.iseg, bp.jseg, 0, &cd2_data.midpoint, &cd2_data.direction);
        for k in 0..3 {
            for j in 0..3 {
                y_full[[k, j]] = y_ref[[k, j]];
            }
        }
    }
    for j in 0..prot2.nres {
        for k in 0..3 {
            y_full[[k, j + 3]] = prot2.ca[[k, j]];
        }
    }
    twist(&mut y_full, 3 + prot2.nres);

    let mut yca = Array2::zeros((3, prot2.nres));
    for j in 0..prot2.nres {
        for k in 0..3 {
            yca[[k, j]] = y_full[[k, j + 3]];
        }
    }

    // Iterative alignment refinement
    let fitz_result = fitz(&mut xca, &yca, RCUT_FITZ, MAXITER);

    // Convert alignment to blocks (1-based residue indices)
    let mut blocks: Vec<AlignmentBlock> = Vec::new();
    for i in 0..prot1.nres {
        if fitz_result.ali[i] >= 0 {
            let j = fitz_result.ali[i] as u32 + 1; // convert to 1-based
            blocks.push(AlignmentBlock {
                l1: i as u32 + 1,
                r1: i as u32 + 1,
                l2: j,
                r2: j,
            });
        }
    }

    if blocks.is_empty() {
        return None;
    }

    // Compress blocks
    let blocks = compress_blocks(&blocks);

    Some(WolfResult {
        cd1: prot1.code.clone(),
        cd2: prot2.code.clone(),
        blocks,
    })
}

/// Run WOLF comparison for all pairs of structures.
///
/// Mimics the Fortran serialcompare loop: for each cd1, build grid once,
/// then compare against all cd2s.
pub fn run_wolf_all_pairs<P: AsRef<Path>>(
    structures: &[&str],
    dat_dir: P,
) -> Vec<WolfResult> {
    let dat_dir = dat_dir.as_ref();
    let mut results = Vec::new();

    // Pre-load all protein data
    let mut protein_data: Vec<Option<WolfProteinData>> = Vec::new();
    for code in structures {
        let filepath = dat_dir.join(format!("{}.dat", code));
        match setup_protein(&filepath) {
            Ok(data) => protein_data.push(Some(data)),
            Err(_) => protein_data.push(None),
        }
    }

    for (_i, cd1_data_opt) in protein_data.iter().enumerate() {
        let cd1_data = match cd1_data_opt {
            Some(d) => d,
            None => continue,
        };

        if cd1_data.protein.nseg <= 2 {
            continue;
        }

        // Build grid for cd1
        let mut grid = SpatialHashGrid::new();
        load_protein(&mut grid, cd1_data);

        for cd2_data_opt in protein_data.iter() {
            let cd2_data = match cd2_data_opt {
                Some(d) => d,
                None => continue,
            };

            if let Some(result) = wolf_compare(cd1_data, cd2_data, &grid) {
                results.push(result);
            }
        }
    }

    results
}
