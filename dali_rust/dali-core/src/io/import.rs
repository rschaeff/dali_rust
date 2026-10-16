//! Full import pipeline: PDB → Protein, as DaliLite's import.pl builds a .dat.
//!
//! PDB-format input follows mpidali.pm prepare_brkfile step for step:
//!   1. dsspcmbi (io::dsspcmbi)            — residues, chain breaks, SS, bridge partners
//!   2. puu on getcoor atoms + that DSSP    — unfolding tree (io::puu)
//!   3. puutos getsecstr                    — SSE segments (io::secstr)
//!   4. puutos prunetree                    — .dat domain trees (io::puutos)
//! CA coordinates come from the DSSP record at 0.1 A, which is all the Fortran
//! pipeline ever sees (.dat stores f8.1).
//!
//! mmCIF input (which DaliLite cannot read) takes the older approximate path:
//! pdbtbx backbone, Rust Kabsch-Sander DSSP, simplified contact-based domain tree.

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
    Dssp(String),
}

impl std::fmt::Display for ImportError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ImportError::Pdb(e) => write!(f, "PDB read error: {}", e),
            ImportError::TooFewResidues(n) => write!(f, "Too few residues: {}", n),
            ImportError::Dssp(e) => write!(f, "DSSP: {}", e),
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
    let name = path.as_ref().to_string_lossy().to_lowercase();
    if name.ends_with(".cif") || name.ends_with(".cif.gz") || name.ends_with(".mmcif") {
        return import_pdb_approximate(path, chain, pdb_code);
    }
    import_pdb_dalilite(path, chain, pdb_code)
}

/// import.pl-equivalent import of a PDB-format file.
pub fn import_pdb_dalilite<P: AsRef<Path>>(
    path: P,
    chain: &str,
    pdb_code: &str,
) -> Result<Protein, ImportError> {
    use super::dsspcmbi::dssp_from_pdb;
    use super::puu::{puu, read_pdb_atoms, to_domain_tree, PuuParams};

    let dssp = dssp_from_pdb(path.as_ref()).map_err(|e| ImportError::Dssp(e.to_string()))?;
    // chain: as given, else the first residue's
    let chainid: u8 = match chain.as_bytes().first() {
        Some(&c) => c,
        None => dssp.res.iter().skip(1).find(|r| !r.is_break()).map(|r| r.chain_id()).unwrap_or(b' '),
    };
    let rows: Vec<&super::dsspcmbi::DsspResidue> = dssp.res.iter().skip(1)
        .filter(|r| !r.is_break() && r.chain_id() == chainid)
        .collect();
    let nres = rows.len();
    if nres < 2 {
        return Err(ImportError::TooFewResidues(nres));
    }
    let dssp_in = dssp.to_puu_input(chainid);

    let atoms = read_pdb_atoms(path.as_ref(), chainid, true)
        .map_err(|e| ImportError::Pdb(PdbReadError::Parse(e.to_string())))?;
    let tree = puu(&atoms, &dssp_in, &PuuParams::default());

    let (segments, secstr, na, nb) = assign_segments(&dssp_in.struc);
    let (domain_tree, _sse_tree) = super::dat::puutos(nres, &segments, &to_domain_tree(&tree));

    let round1 = |x: f64| -> f64 { format!("{:.1}", x).parse().unwrap() };
    let mut ca = Array2::<f64>::zeros((3, nres));
    let mut sequence = String::with_capacity(nres);
    let mut resid_map = Vec::with_capacity(nres);
    for (j, r) in rows.iter().enumerate() {
        for k in 0..3 {
            ca[[k, j]] = round1(r.ca[k]);
        }
        sequence.push(r.aa as char);
        resid_map.push(String::from_utf8_lossy(&r.aaident[0..4]).trim().parse().unwrap_or(0));
    }

    let chain_str = (chainid as char).to_string();
    let code = if pdb_code.len() <= 4 && chainid != b' ' && !pdb_code.ends_with(&*chain_str) {
        format!("{}{}", pdb_code, chain_str)
    } else {
        pdb_code.to_string()
    };

    Ok(Protein {
        code,
        nres,
        nseg: segments.len(),
        na,
        nb,
        segments,
        secstr,
        ca,
        sequence,
        domain_tree,
        resid_map,
        numbering: crate::ResidNumbering::Pdb,
    })
}

/// Approximate import (mmCIF): pdbtbx backbone, Rust DSSP, simplified domain tree.
pub fn import_pdb_approximate<P: AsRef<Path>>(
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
