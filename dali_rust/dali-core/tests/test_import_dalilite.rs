//! import_pdb from a PDB file alone reproduces DaliLite import.pl's .dat
//! (sections 1-3: header, SSE segments, CA, SSE tree, residue tree) byte for byte.

use dali_core::io::{import_pdb, write_dat};

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

/// Lines of sections 1-3 (everything before the first '-' record).
fn sections(text: &str) -> Vec<String> {
    text.lines().take_while(|l| !l.starts_with('-')).map(|l| l.trim_end().to_string()).collect()
}

#[test]
fn test_import_reproduces_fortran_dat() {
    let mut ids: Vec<String> = std::fs::read_dir(format!("{}/puu", V)).unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.path().join("units.puu").exists())
        .map(|e| e.file_name().to_string_lossy().to_string())
        .collect();
    ids.sort();
    assert!(ids.len() >= 28);
    let out = std::env::temp_dir().join(format!("dali_import_{}", std::process::id()));
    std::fs::create_dir_all(&out).unwrap();
    let mut bad = Vec::new();
    for id in &ids {
        let p = import_pdb(pdb_path(id), &id[4..5], &id[..4]).unwrap();
        assert_eq!(p.code, *id);
        let dst = out.join(format!("{}.dat", id));
        write_dat(&p, &dst).unwrap();
        let ours = sections(&std::fs::read_to_string(&dst).unwrap());
        let theirs = sections(&std::fs::read_to_string(format!("{}/puu/{}/{}.dat", V, id, id)).unwrap());
        if ours != theirs {
            let k = ours.iter().zip(&theirs).position(|(a, b)| a != b).unwrap_or(ours.len().min(theirs.len()));
            eprintln!("{}: {} vs {} lines, first diff line {}:\n  ours    |{}|\n  fortran |{}|", id, ours.len(), theirs.len(), k + 1,
                ours.get(k).map(|s| s.as_str()).unwrap_or("-"), theirs.get(k).map(|s| s.as_str()).unwrap_or("-"));
            bad.push(id.clone());
        }
    }
    let _ = std::fs::remove_dir_all(&out);
    eprintln!("import: {}/{} .dat byte-identical (sections 1-3)", ids.len() - bad.len(), ids.len());
    assert!(bad.is_empty(), "{:?}", bad);
}
