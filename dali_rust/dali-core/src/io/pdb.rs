use std::io::{BufReader, Read as IoRead};
use std::path::Path;

/// Per-residue backbone atom coordinates.
#[derive(Debug, Clone)]
pub struct BackboneResidue {
    pub n: [f64; 3],
    pub ca: [f64; 3],
    pub c: [f64; 3],
    pub o: [f64; 3],
    pub resname: String, // 3-letter code e.g. "ALA"
}

/// Backbone data extracted from a PDB/CIF file.
#[derive(Debug, Clone)]
pub struct BackboneData {
    pub code: String,
    pub chain_id: String,
    pub residues: Vec<BackboneResidue>,
    pub sequence: String, // 1-letter codes
}

/// Errors from PDB reading.
#[derive(Debug)]
pub enum PdbReadError {
    Io(String),
    Parse(String),
    NoChain(String),
}

impl std::fmt::Display for PdbReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PdbReadError::Io(msg) => write!(f, "PDB I/O error: {}", msg),
            PdbReadError::Parse(msg) => write!(f, "PDB parse error: {}", msg),
            PdbReadError::NoChain(msg) => write!(f, "PDB chain not found: {}", msg),
        }
    }
}

impl std::error::Error for PdbReadError {}

/// Convert 3-letter amino acid code to 1-letter code.
fn aa3to1(resname: &str) -> char {
    match resname {
        "ALA" => 'A', "ARG" => 'R', "ASN" => 'N', "ASP" => 'D',
        "CYS" => 'C', "GLN" => 'Q', "GLU" => 'E', "GLY" => 'G',
        "HIS" => 'H', "ILE" => 'I', "LEU" => 'L', "LYS" => 'K',
        "MET" => 'M', "PHE" => 'F', "PRO" => 'P', "SER" => 'S',
        "THR" => 'T', "TRP" => 'W', "TYR" => 'Y', "VAL" => 'V',
        "MSE" => 'M', // selenomethionine
        _ => 'X',
    }
}

/// Detect file format from path, ignoring .gz extension.
fn detect_format(path: &str) -> pdbtbx::Format {
    let lower = path.to_lowercase();
    let base = if lower.ends_with(".gz") {
        &lower[..lower.len() - 3]
    } else {
        &lower
    };
    if base.ends_with(".cif") || base.ends_with(".mmcif") {
        pdbtbx::Format::Mmcif
    } else {
        // Default to PDB for .pdb, .ent, or anything else
        pdbtbx::Format::Pdb
    }
}

/// Read a PDB or CIF file and extract backbone data for the specified chain.
///
/// - `path`: Path to PDB/CIF file (supports .gz compression, including .ent.gz)
/// - `chain`: Chain identifier (e.g. "A"). If empty, uses the first chain.
/// - `pdb_code`: 4-character PDB code (e.g. "101m")
///
/// Returns `BackboneData` with per-residue N, CA, C, O coordinates and sequence.
pub fn read_pdb<P: AsRef<Path>>(
    path: P,
    chain: &str,
    pdb_code: &str,
) -> Result<BackboneData, PdbReadError> {
    let path_ref = path.as_ref();
    let path_str = path_ref.to_str().unwrap_or("");
    let is_gz = path_str.ends_with(".gz");
    let format = detect_format(path_str);

    // Read and decompress file content
    let file = std::fs::File::open(path_ref)
        .map_err(|e| PdbReadError::Io(format!("Cannot open {}: {}", path_str, e)))?;

    let reader: Box<dyn IoRead> = if is_gz {
        Box::new(flate2::read::GzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let buf_reader = BufReader::new(reader);

    let mut opts = pdbtbx::ReadOptions::new();
    opts.set_format(format)
        .set_only_first_model(true)
        .set_only_atomic_coords(true)
        .set_level(pdbtbx::StrictnessLevel::Loose);

    let (pdb, _errors) = opts.read_raw(buf_reader).map_err(|errs| {
        PdbReadError::Parse(format!(
            "Failed to parse {}: {:?}",
            path_str,
            errs.iter().map(|e| format!("{}", e)).collect::<Vec<_>>()
        ))
    })?;

    extract_backbone(&pdb, chain, pdb_code)
}

/// Extract backbone data from a parsed PDB structure.
fn extract_backbone(
    pdb: &pdbtbx::PDB,
    chain: &str,
    pdb_code: &str,
) -> Result<BackboneData, PdbReadError> {
    // Find the target chain
    let target_chain = if chain.is_empty() {
        pdb.chains().next().ok_or_else(|| {
            PdbReadError::NoChain("No chains found in PDB".to_string())
        })?
    } else {
        pdb.chains()
            .find(|c| c.id() == chain)
            .ok_or_else(|| {
                PdbReadError::NoChain(format!("Chain '{}' not found", chain))
            })?
    };

    let chain_id = target_chain.id().to_string();

    let mut residues = Vec::new();
    let mut sequence = String::new();

    for residue in target_chain.residues() {
        let resname = match residue.name() {
            Some(n) => n.to_string(),
            None => continue,
        };

        // Find backbone atoms N, CA, C, O
        let mut n_coord: Option<[f64; 3]> = None;
        let mut ca_coord: Option<[f64; 3]> = None;
        let mut c_coord: Option<[f64; 3]> = None;
        let mut o_coord: Option<[f64; 3]> = None;

        for atom in residue.atoms() {
            let name = atom.name().trim().to_uppercase();
            let (x, y, z) = atom.pos();
            match name.as_str() {
                "N" => n_coord = Some([x, y, z]),
                "CA" => ca_coord = Some([x, y, z]),
                "C" => c_coord = Some([x, y, z]),
                "O" => o_coord = Some([x, y, z]),
                _ => {}
            }
        }

        // Skip residues with missing backbone atoms
        let (n, ca, c, o) = match (n_coord, ca_coord, c_coord, o_coord) {
            (Some(n), Some(ca), Some(c), Some(o)) => (n, ca, c, o),
            _ => continue,
        };

        // Only include amino acids
        let aa1 = aa3to1(&resname);
        if aa1 == 'X' && resname != "UNK" {
            continue;
        }

        sequence.push(aa1);
        residues.push(BackboneResidue {
            n,
            ca,
            c,
            o,
            resname,
        });
    }

    let code = format!(
        "{}{}",
        pdb_code,
        if chain_id.len() == 1 { &chain_id } else { "" }
    );

    Ok(BackboneData {
        code,
        chain_id,
        residues,
        sequence,
    })
}
