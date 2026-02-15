//! Full import pipeline: PDB/CIF → Protein.
//!
//! Reads a PDB or mmCIF file and produces a `Protein` struct ready for the
//! DALI comparison pipeline. Steps:
//!   1. Parse backbone atoms (N, CA, C, O) via pdbtbx
//!   2. Compute 8-state DSSP secondary structure (Kabsch-Sander)
//!   3. Reduce to 3-state and assign SSE segments (puutos.f port)
//!   4. Compute hierarchical domain decomposition (puu.f port)
//!   5. Assemble into Protein struct

use std::path::Path;
use ndarray::Array2;

use crate::types::Protein;
use super::pdb::{read_pdb, PdbReadError};
use super::dssp::compute_dssp;
use super::secstr::assign_segments;
use super::domain::build_domain_tree;

/// Errors from the import pipeline.
#[derive(Debug)]
pub enum ImportError {
    Pdb(PdbReadError),
    TooFewResidues(usize),
}

impl std::fmt::Display for ImportError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ImportError::Pdb(e) => write!(f, "PDB read error: {}", e),
            ImportError::TooFewResidues(n) => write!(f, "Too few residues: {}", n),
        }
    }
}

impl std::error::Error for ImportError {}

impl From<PdbReadError> for ImportError {
    fn from(e: PdbReadError) -> Self {
        ImportError::Pdb(e)
    }
}

/// Import a PDB/CIF file and produce a Protein struct.
///
/// - `path`: Path to PDB/CIF file (supports .gz compression)
/// - `chain`: Chain identifier (e.g. "A"). If empty, uses the first chain.
/// - `pdb_code`: 4-character PDB code (e.g. "101m")
pub fn import_pdb<P: AsRef<Path>>(
    path: P,
    chain: &str,
    pdb_code: &str,
) -> Result<Protein, ImportError> {
    // Step 1: Read backbone atoms
    let backbone = read_pdb(path, chain, pdb_code)?;
    let nres = backbone.residues.len();

    if nres < 2 {
        return Err(ImportError::TooFewResidues(nres));
    }

    // Step 2: Compute 8-state DSSP
    let dssp8 = compute_dssp(&backbone.residues);

    // Step 3: Reduce to 3-state + assign SSE segments
    let (segments, secstr, na, nb) = assign_segments(&dssp8);
    let nseg = segments.len();

    // Step 4: Build CA coordinate array (3, nres)
    let mut ca = Array2::<f64>::zeros((3, nres));
    for (j, res) in backbone.residues.iter().enumerate() {
        ca[[0, j]] = res.ca[0];
        ca[[1, j]] = res.ca[1];
        ca[[2, j]] = res.ca[2];
    }

    // Step 5: Compute domain decomposition
    let ca_flat: Vec<[f64; 3]> = backbone.residues.iter().map(|r| r.ca).collect();
    let domain_tree = build_domain_tree(&ca_flat, &segments);

    Ok(Protein {
        code: backbone.code,
        nres,
        nseg,
        na,
        nb,
        segments,
        secstr,
        ca,
        sequence: backbone.sequence,
        domain_tree,
        resid_map: backbone.resid_map,
        numbering: crate::ResidNumbering::Pdb,
    })
}
