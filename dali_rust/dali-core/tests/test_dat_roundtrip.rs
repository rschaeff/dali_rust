//! .dat round-trip: what Rust writes must be readable by every consumer, PARSI included.
//!
//! The module tests all read Fortran-written fixtures with 5-char codes, so until
//! 2026-10 nothing exercised PARSI on a Rust-written .dat. It could not read one: the
//! header was whitespace-delimited where PARSI expected fixed columns, and the SSE-level
//! domain tree PARSI needs was never written. Every Rust-imported query was skipped
//! silently. These tests pin the formats against Fortran output.

use std::fs;
use std::path::PathBuf;

use dali_core::io::{derive_sse_tree, import_pdb, read_dat, write_dat};
use dali_core::parsi::protein_io::parsireadproteindata;
use dali_core::pipeline::run_parsi_path;
use dali_core::ProteinStore;

const DAT_DIR_5: &str = concat!(
    env!("CARGO_MANIFEST_DIR"), "/../..", "/validation/fixtures/structures"
);
const PDB_DIR_5: &str = concat!(
    env!("CARGO_MANIFEST_DIR"), "/../..", "/validation/fixtures/structures/pdb"
);
const DAT_DIR_37: &str = concat!(
    env!("CARGO_MANIFEST_DIR"), "/../..", "/validation/corpus_expansion/dat"
);

fn scratch(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("dali_roundtrip_{}_{}", name, std::process::id()));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    d
}

fn fortran_dats() -> Vec<PathBuf> {
    let mut v = Vec::new();
    for dir in [DAT_DIR_5, DAT_DIR_37] {
        for e in fs::read_dir(dir).unwrap() {
            let p = e.unwrap().path();
            if p.extension().map_or(false, |x| x == "dat") {
                v.push(p);
            }
        }
    }
    v.sort();
    v
}

/// Lines of the n-th `>>>>` section (1-based), header included.
fn section(text: &str, n: usize) -> Vec<String> {
    let mut k = 0;
    let mut out = Vec::new();
    for line in text.lines() {
        if line.starts_with(">>>>") {
            k += 1;
        }
        if k == n {
            out.push(line.to_string());
        } else if k > n {
            break;
        }
    }
    out
}

/// Section 2 as Fortran wrote it: (child1, child2, seglist) per node.
fn fortran_sse_tree(text: &str) -> Vec<(usize, usize, Vec<usize>)> {
    section(text, 2)[1..]
        .iter()
        .map(|l| {
            let f = |a: usize, b: usize| l[a..b].trim().parse::<usize>().unwrap();
            let dns = f(15, 19);
            let segs = (0..dns).map(|j| f(19 + 4 * j, 23 + 4 * j)).collect();
            (f(7, 11), f(11, 15), segs)
        })
        .collect()
}

/// derive_sse_tree applied to Fortran's residue-level tree (section 3) must reproduce
/// Fortran's SSE-level tree (section 2): same puutos.f logic, same input.
#[test]
fn test_derive_sse_tree_matches_fortran() {
    let files = fortran_dats();
    assert!(files.len() >= 42, "expected the 5 + 37 Fortran fixtures, found {}", files.len());
    let mut bad = Vec::new();
    for path in &files {
        let text = fs::read_to_string(path).unwrap();
        let p = read_dat(path).unwrap();
        let derived: Vec<(usize, usize, Vec<usize>)> =
            derive_sse_tree(p.nres, &p.segments, &p.domain_tree)
                .into_iter()
                .map(|n| (n.left_child, n.right_child, n.sses))
                .collect();
        if derived != fortran_sse_tree(&text) {
            bad.push(path.file_name().unwrap().to_string_lossy().to_string());
        }
    }
    assert!(bad.is_empty(), "SSE tree differs from Fortran for {}/{}: {:?}",
            bad.len(), files.len(), bad);
}

/// For a 5-char code, write_dat(read_dat(f)) reproduces Fortran sections 1 and 2
/// byte for byte, apart from the node-type column.
#[test]
fn test_write_dat_reproduces_fortran_sections_1_2() {
    let out = scratch("sections");
    let mut bad = Vec::new();
    for path in fortran_dats() {
        let text = fs::read_to_string(&path).unwrap();
        let p = read_dat(&path).unwrap();
        let dst = out.join(path.file_name().unwrap());
        write_dat(&p, &dst).unwrap();
        let ours = fs::read_to_string(&dst).unwrap();
        // Column 5 is the node type; NodeType folds Fortran's '=', 'e', 'n' into '-'.
        // No consumer reads it from section 2 (DP reads '+'/'*' from section 3).
        let mask = |v: Vec<String>| -> Vec<String> {
            v.into_iter()
                .map(|l| if l.starts_with(">>>>") || l.len() < 6 { l }
                         else { format!("{}?{}", &l[..5], &l[6..]) })
                .collect()
        };
        for n in [1, 2] {
            if mask(section(&ours, n)) != mask(section(&text, n)) {
                bad.push(format!("{} section {}", path.file_name().unwrap().to_string_lossy(), n));
            }
        }
    }
    let _ = fs::remove_dir_all(&out);
    assert!(bad.is_empty(), "{} sections differ: {:?}", bad.len(), &bad[..bad.len().min(10)]);
}

/// Section 3 survives a write/read cycle (children, segments) for Fortran input.
#[test]
fn test_residue_tree_roundtrip() {
    let out = scratch("restree");
    for path in fortran_dats() {
        let p = read_dat(&path).unwrap();
        let dst = out.join(path.file_name().unwrap());
        write_dat(&p, &dst).unwrap();
        let q = read_dat(&dst).unwrap();
        assert_eq!(p.domain_tree.len(), q.domain_tree.len());
        for (a, b) in p.domain_tree.iter().zip(&q.domain_tree) {
            assert_eq!((a.left_child, a.right_child, &a.segments),
                       (b.left_child, b.right_child, &b.segments),
                       "{:?} node {}", path, a.index);
        }
    }
    let _ = fs::remove_dir_all(&out);
}

/// Fortran section-3 child pointers sit at columns 7-14; they used to be read from 6-18
/// in width-4 chunks, which garbled them.
#[test]
fn test_read_dat_children_fortran_columns() {
    let p = read_dat(format!("{}/1a87A.dat", DAT_DIR_5)).unwrap();
    let text = fs::read_to_string(format!("{}/1a87A.dat", DAT_DIR_5)).unwrap();
    for (node, line) in p.domain_tree.iter().zip(&section(&text, 3)[1..]) {
        let c1: usize = line[7..11].trim().parse().unwrap();
        let c2: usize = line[11..15].trim().parse().unwrap();
        assert_eq!((node.left_child, node.right_child), (c1, c2), "{}", line);
    }
}

/// The DPAM shape: Rust-imported proteins with codes longer than 5 characters. PARSI
/// must read them and find the 101m/1a00 globin match (Z 20.6 from Fortran .dat).
#[test]
fn test_parsi_on_rust_imported_long_codes() {
    let dir = scratch("longcode");
    let store = ProteinStore::new(&dir);
    let q = import_pdb(format!("{}/pdb101m.ent.gz", PDB_DIR_5), "A", "query_101m").unwrap();
    let t = import_pdb(format!("{}/pdb1a00.ent.gz", PDB_DIR_5), "A", "000140533").unwrap();
    let (qc, tc) = (q.code.clone(), t.code.clone());
    assert!(qc.len() > 5 && tc.len() > 5, "codes {} {}", qc, tc);
    store.add_protein(q).unwrap();
    store.add_protein(t).unwrap();

    for c in [&qc, &tc] {
        let path = format!("{}/{}.dat", dir.display(), c);
        let pp = parsireadproteindata(&path).unwrap_or_else(|e| panic!("{}: {}", path, e));
        assert_eq!(&pp.code, c);
        assert!(pp.ndom >= 1);
    }

    for (a, b) in [(&qc, &tc), (&tc, &qc)] {
        let hits = run_parsi_path(a, &[b.as_str()], &store);
        let best = hits.iter().map(|h| h.zscore).fold(f64::MIN, f64::max);
        assert!(best > 15.0, "PARSI {} -> {}: best Z {} from {} hits", a, b, best, hits.len());
    }
    let _ = fs::remove_dir_all(&dir);
}

/// Files from the pre-2026-10 writer (whitespace header, empty section 2) are still
/// readable: PARSI derives the SSE tree from section 3.
#[test]
fn test_parsi_reads_legacy_rust_dat() {
    let dir = scratch("legacy");
    let p = read_dat(format!("{}/1a87A.dat", DAT_DIR_5)).unwrap();
    let new = dir.join("new.dat");
    write_dat(&p, &new).unwrap();
    let text = fs::read_to_string(&new).unwrap();

    // Old layout: whitespace headers, section 2 = bare copy of the section 1 header.
    let secstr: String = p.secstr.iter().map(|s| s.to_char()).collect();
    let h1 = format!(">>>> {} {} {} {} {} {}", p.code, p.nres, p.nseg, p.na, p.nb, secstr);
    let mut legacy = vec![h1.clone()];
    legacy.extend(section(&text, 1).into_iter().skip(1));
    legacy.push(h1);
    let s3 = section(&text, 3);
    legacy.push(format!(">>>> {} {}", p.code, p.domain_tree.len()));
    legacy.extend(s3.into_iter().skip(1));
    let old = dir.join("old.dat");
    fs::write(&old, legacy.join("\n") + "\n").unwrap();

    let a = parsireadproteindata(new.to_str().unwrap()).unwrap();
    let b = parsireadproteindata(old.to_str().unwrap()).unwrap();
    assert_eq!((a.ndom, &a.node_child, &a.domns, &a.domseglist),
               (b.ndom, &b.node_child, &b.domns, &b.domseglist));
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn test_parsi_read_error_is_reported() {
    let dir = scratch("err");
    let bad = dir.join("bad.dat");
    fs::write(&bad, ">>>> x notanumber\n").unwrap();
    assert!(parsireadproteindata(bad.to_str().unwrap()).is_err());
    assert!(parsireadproteindata(dir.join("missing.dat").to_str().unwrap()).is_err());
    let _ = fs::remove_dir_all(&dir);
}

/// Fortran singletex treats a segment longer than 100 residues as fatal (its score
/// table is 100x100): lskipflag is set and PARSI abandons the rest of the pair,
/// keeping hits already written. AlphaFold models often have such helices. The port
/// used to return silently and keep searching with no end-trimming for that segment,
/// emitting alignments past the end of cd2 that crashed FILTER95.
#[test]
fn test_parsi_long_segment_aborts_like_fortran() {
    use dali_core::parsi::run_parsi;
    use dali_core::{Segment, SseType};

    let dir = scratch("longseg");
    let mut p = read_dat(format!("{}/1a87A.dat", DAT_DIR_5)).unwrap();
    // Replace everything up to residue 110 with one 110-residue helix.
    let mut segs = vec![Segment {
        start: 1, end: 110, sse_type: SseType::Helix, check_start: 1, check_end: 110, checkx: 0,
    }];
    segs.extend(p.segments.iter().filter(|s| s.start > 110).cloned());
    p.secstr = segs.iter().map(|s| s.sse_type).collect();
    p.na = p.secstr.iter().filter(|&&s| s == SseType::Helix).count();
    p.nb = p.secstr.len() - p.na;
    p.nseg = segs.len();
    p.segments = segs;
    p.code = "longseg".into();
    write_dat(&p, dir.join("longseg.dat")).unwrap();
    fs::copy(format!("{}/101mA.dat", DAT_DIR_5), dir.join("101mA.dat")).unwrap();

    let tree = derive_sse_tree(p.nres, &p.segments, &p.domain_tree);
    let hits = run_parsi(&["longseg".to_string(), "101mA".to_string()], dir.to_str().unwrap(), None);
    let mut n = 0;
    for h in hits.iter().filter(|h| h.cd1 == "longseg" && h.cd2 == "101mA") {
        n += 1;
        assert!(h.idom > tree.len() || !tree[h.idom - 1].sses.contains(&1),
                "idom {} contains the 110-residue segment", h.idom);
        for (a, b) in h.ranges_cd1.iter().zip(&h.ranges_cd2) {
            if a.0 > 0 {
                assert!(a.1 as usize <= p.nres && b.0 >= 1 && b.1 <= 154, "{:?} -> {:?}", a, b);
            }
        }
    }
    // the full pipeline must not panic on it either
    let store = ProteinStore::new(&dir);
    let _ = run_parsi_path("longseg", &["101mA"], &store);
    eprintln!("longseg -> 101mA: {} PARSI hits from domains without the long segment", n);
    let _ = fs::remove_dir_all(&dir);
}

/// PARSI results must not depend on code length. Codes shorter than 5 characters
/// used to be space-padded into the refine-line tag, which the parser then split,
/// silently dropping every hit -- DPAM's RustDALI names its query "q" -> "qA".
#[test]
fn test_parsi_independent_of_code_length() {
    let dir = scratch("codelen");
    let mut results = Vec::new();
    for (a, b) in [("qA", "tA"), ("101mA", "1a00A"), ("query_101m", "000140533")] {
        fs::copy(format!("{}/101mA.dat", DAT_DIR_5), dir.join(format!("{}.dat", a))).unwrap();
        fs::copy(format!("{}/1a00A.dat", DAT_DIR_5), dir.join(format!("{}.dat", b))).unwrap();
        let store = ProteinStore::new(&dir);
        let mut z: Vec<i64> = run_parsi_path(a, &[b], &store).iter()
            .map(|h| (h.zscore * 1000.0).round() as i64).collect();
        z.sort();
        assert!(!z.is_empty(), "{} -> {}: no PARSI hits", a, b);
        results.push(z);
    }
    assert_eq!(results[0], results[1]);
    assert_eq!(results[1], results[2]);
    let _ = fs::remove_dir_all(&dir);
}
