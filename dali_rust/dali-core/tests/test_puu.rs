//! puu port vs Fortran puu, fed the same PDB and the same (Fortran dsspcmbi) DSSP,
//! so DSSP divergence is excluded. Reference: validation/puu/<id>/units.puu, made by
//! validation/puu/make_puu_fixtures.sh with DaliLite's own binaries.

use dali_core::io::puu::{puu, read_dssp, read_pdb_atoms, read_units_puu, PuuNode, PuuParams};

const V: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation");

fn pdb_path(id: &str) -> String {
    let code = &id[..4];
    let five = ["101mA", "1a00A", "1a87A", "1allA", "1binA"];
    if five.contains(&id) {
        format!("{}/fixtures/structures/pdb/pdb{}.ent.gz", V, code)
    } else {
        format!("{}/corpus_expansion/pdb/pdb{}.ent.gz", V, code)
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

fn same(a: &PuuNode, b: &PuuNode) -> bool {
    let r3 = |x: f32| (x.min(9999.999) * 1000.0).round();
    a.child == b.child && a.node_type == b.node_type && a.size == b.size
        && a.nbp == b.nbp && a.ranges == b.ranges
        && r3(a.tau) == r3(b.tau) && r3(a.tert) == r3(b.tert)
}

#[test]
fn test_puu_matches_fortran() {
    let ids = ids();
    assert!(ids.len() >= 23, "{} fixtures", ids.len());
    let mut bad = Vec::new();
    for id in &ids {
        let chain = id.as_bytes()[4];
        let dir = format!("{}/puu/{}", V, id);
        let s = match std::fs::read_dir(&dir).unwrap().filter_map(|e| e.ok())
            .find(|e| e.file_name().to_string_lossy().ends_with(".pdb.gz")) {
            Some(e) => read_pdb_atoms(e.path(), chain, true).unwrap(),
            None => read_pdb_atoms(pdb_path(id), chain, true).unwrap(),
        };
        let dssp_file = std::fs::read_dir(&dir).unwrap().filter_map(|e| e.ok())
            .find(|e| e.file_name().to_string_lossy().ends_with(".dssp")).unwrap().path();
        let d = read_dssp(&dssp_file, chain).unwrap();
        let ours = puu(&s, &d, &PuuParams::default());
        let theirs = read_units_puu(format!("{}/units.puu", dir)).unwrap();
        let first_diff = ours.iter().zip(&theirs).position(|(a, b)| !same(a, b));
        if ours.len() != theirs.len() || first_diff.is_some() {
            let k = first_diff.unwrap_or(ours.len().min(theirs.len()));
            eprintln!("{}: nres {} nodes ours {} fortran {}; first diff node {}:\n  ours    {:?}\n  fortran {:?}",
                id, s.nres, ours.len(), theirs.len(), k + 1, ours.get(k), theirs.get(k));
            bad.push(id.clone());
        }
    }
    eprintln!("puu: {}/{} identical", ids.len() - bad.len(), ids.len());
    assert!(bad.is_empty(), "differs: {:?}", bad);
}

/// puu -> getsecstr -> puutos reproduces Fortran's .dat segments (section 1), SSE tree
/// (section 2) and pruned residue tree (section 3), from the same PDB and DSSP.
#[test]
fn test_puu_puutos_matches_fortran_dat() {
    use dali_core::io::{assign_segments, puutos, read_dat};
    use dali_core::io::puu::to_domain_tree;

    let mut bad = Vec::new();
    let ids = ids();
    for id in &ids {
        let chain = id.as_bytes()[4];
        let dir = format!("{}/puu/{}", V, id);
        let pdb = std::fs::read_dir(&dir).unwrap().filter_map(|e| e.ok())
            .find(|e| e.file_name().to_string_lossy().ends_with(".pdb.gz"))
            .map(|e| e.path().to_string_lossy().to_string())
            .unwrap_or_else(|| pdb_path(id));
        let s = read_pdb_atoms(&pdb, chain, true).unwrap();
        let dssp_file = std::fs::read_dir(&dir).unwrap().filter_map(|e| e.ok())
            .find(|e| e.file_name().to_string_lossy().ends_with(".dssp")).unwrap().path();
        let d = read_dssp(&dssp_file, chain).unwrap();
        let tree = puu(&s, &d, &PuuParams::default());

        let (segments, _secstr, _na, _nb) = assign_segments(&d.struc);
        let (res_tree, sse_tree) = puutos(s.nres, &segments, &to_domain_tree(&tree));
        let f = read_dat(format!("{}/{}.dat", dir, id)).unwrap();
        let ftext = std::fs::read_to_string(format!("{}/{}.dat", dir, id)).unwrap();

        let segs_ok = segments.len() == f.segments.len()
            && segments.iter().zip(&f.segments).all(|(a, b)| (a.start, a.end, a.check_start, a.check_end, a.checkx)
                == (b.start, b.end, b.check_start, b.check_end, b.checkx));
        let res_ok = res_tree.len() == f.domain_tree.len()
            && res_tree.iter().zip(&f.domain_tree).all(|(a, b)|
                (a.left_child, a.right_child, &a.segments, a.node_type) == (b.left_child, b.right_child, &b.segments, b.node_type));
        // section 2 as Fortran wrote it
        let mut k = 0;
        let fsse: Vec<(usize, usize, Vec<usize>)> = ftext.lines()
            .filter(|l| { if l.starts_with(">>>>") { k += 1; false } else { k == 2 } })
            .map(|l| {
                let g = |a: usize, b: usize| l[a..b].trim().parse::<usize>().unwrap();
                let n = g(15, 19);
                (g(7, 11), g(11, 15), (0..n).map(|j| g(19 + 4 * j, 23 + 4 * j)).collect())
            })
            .collect();
        let ours: Vec<(usize, usize, Vec<usize>)> = sse_tree.iter()
            .map(|n| (n.left_child, n.right_child, n.sses.clone())).collect();
        let sse_ok = ours == fsse;
        if !(segs_ok && res_ok && sse_ok) {
            eprintln!("{}: segments {} section3 {} section2 {}", id, segs_ok, res_ok, sse_ok);
            bad.push(id.clone());
        }
    }
    eprintln!("puu+puutos: {}/{} match Fortran .dat", ids.len() - bad.len(), ids.len());
    assert!(bad.is_empty(), "{:?}", bad);
}
