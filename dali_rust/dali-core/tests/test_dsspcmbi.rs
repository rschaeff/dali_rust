//! dsspcmbi port vs DaliLite's dsspcmbi binary, line by line, on the puu fixtures
//! (validation/puu/<id>/<code>.dssp were written by the Fortran pipeline's own
//! dsspcmbi from the same PDB files).

use dali_core::io::dsspcmbi::{dssp_from_pdb, Dssp};

const V: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation");

fn pdb_path(id: &str) -> String {
    let dir = format!("{}/puu/{}", V, id);
    if let Some(e) = std::fs::read_dir(&dir).unwrap().filter_map(|e| e.ok())
        .find(|e| e.file_name().to_string_lossy().ends_with(".pdb.gz")) {
        return e.path().to_string_lossy().to_string();
    }
    let five = ["101mA", "1a00A", "1a87A", "1allA", "1binA"];
    if five.contains(&id) {
        format!("{}/fixtures/structures/pdb/pdb{}.ent.gz", V, &id[..4])
    } else {
        format!("{}/corpus_expansion/pdb/pdb{}.ent.gz", V, &id[..4])
    }
}

fn ids() -> Vec<String> {
    let mut v: Vec<String> = std::fs::read_dir(format!("{}/puu", V)).unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.path().join("units.puu").exists())
        .map(|e| e.file_name().to_string_lossy().to_string())
        .collect();
    v.sort();
    v
}

/// Our residue i as the DSSP columns we port: 1-33 minus accessibility and H-bond
/// columns: "#####_RRRRIC_A?_S_TTTBCLL BP1 BP2L"
fn our_line(d: &Dssp, i: usize) -> String {
    let r = &d.res[i];
    let mut s = format!("{:5} {} {}  {} ", i, String::from_utf8_lossy(&r.aaident), r.aa as char, r.ss[0] as char);
    for k in 1..8 {
        s.push(r.ss[k] as char);
    }
    s.push_str(&format!("{:4}{:4}{}", r.partner[0], r.partner[1], r.sheetlabel as char));
    s
}

/// Fortran line normalised the same way (pdb-chain-break column 15 dropped,
/// SS-bridge lowercase cysteine labels mapped back to C).
fn their_line(l: &str) -> String {
    let b = l.as_bytes();
    let mut aa = b[13];
    if aa.is_ascii_lowercase() {
        aa = b'C';
    }
    format!("{} {}  {}", &l[0..12], aa as char, &l[16..34])
}

#[test]
fn test_dsspcmbi_matches_fortran() {
    let ids = ids();
    assert!(ids.len() >= 28);
    let mut bad = Vec::new();
    let mut nres_total = 0;
    for id in &ids {
        let dir = format!("{}/puu/{}", V, id);
        let dssp_file = std::fs::read_dir(&dir).unwrap().filter_map(|e| e.ok())
            .find(|e| e.file_name().to_string_lossy().ends_with(".dssp")).unwrap().path();
        let text = std::fs::read_to_string(&dssp_file).unwrap();
        let theirs: Vec<String> = text.lines()
            .skip_while(|l| !l.starts_with("  #  RESIDUE"))
            .skip(1)
            .map(their_line)
            .collect();
        let d = dssp_from_pdb(pdb_path(id)).unwrap();
        let ours: Vec<String> = (1..d.res.len()).map(|i| our_line(&d, i)).collect();
        nres_total += theirs.len();
        let first = ours.iter().zip(&theirs).position(|(a, b)| a != b);
        if ours.len() != theirs.len() || first.is_some() {
            let k = first.unwrap_or(ours.len().min(theirs.len()));
            let ndiff = ours.iter().zip(&theirs).filter(|(a, b)| a != b).count();
            eprintln!("{}: lines ours {} fortran {}, {} differ; first at {}:\n  ours    |{}|\n  fortran |{}|",
                id, ours.len(), theirs.len(), ndiff, k + 1,
                ours.get(k).map(|s| s.as_str()).unwrap_or("-"), theirs.get(k).map(|s| s.as_str()).unwrap_or("-"));
            bad.push(id.clone());
        }
    }
    eprintln!("dsspcmbi: {}/{} files identical ({} residue lines)", ids.len() - bad.len(), ids.len(), nres_total);
    assert!(bad.is_empty(), "{:?}", bad);
}
