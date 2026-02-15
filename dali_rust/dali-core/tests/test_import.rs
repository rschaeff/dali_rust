use dali_core::io::{read_dat, write_dat, read_pdb, import_pdb};

const DAT_DIR_5: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../..",
    "/validation/fixtures/structures"
);
const PDB_DIR_5: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../..",
    "/validation/fixtures/structures/pdb"
);
const DAT_DIR_18: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../..",
    "/validation/corpus_expansion/ground_truth_18/structures"
);
const PDB_DIR_18: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../..",
    "/validation/corpus_expansion/pdb"
);

/// Test that PDB reader extracts correct nres and CA coordinates for the 5-struct corpus.
#[test]
fn test_pdb_reader_5struct_nres_and_ca() {
    let cases = [
        ("101m", "A", "101mA"),
        ("1a00", "A", "1a00A"),
        ("1a87", "A", "1a87A"),
        ("1all", "A", "1allA"),
        ("1bin", "A", "1binA"),
    ];

    for (pdb_code, chain, dat_code) in &cases {
        let pdb_path = format!("{}/pdb{}.ent.gz", PDB_DIR_5, pdb_code);
        let dat_path = format!("{}/{}.dat", DAT_DIR_5, dat_code);

        let backbone = read_pdb(&pdb_path, chain, pdb_code).unwrap_or_else(|e| {
            panic!("Failed to read PDB {}: {}", pdb_code, e);
        });
        let protein = read_dat(&dat_path).unwrap_or_else(|e| {
            panic!("Failed to read DAT {}: {}", dat_code, e);
        });

        assert_eq!(
            backbone.residues.len(),
            protein.nres,
            "{}: nres mismatch (PDB={}, DAT={})",
            dat_code,
            backbone.residues.len(),
            protein.nres
        );

        // Verify CA coordinates match within tolerance
        let tol = 0.15;
        for (i, res) in backbone.residues.iter().enumerate() {
            let dx = (res.ca[0] - protein.ca[[0, i]]).abs();
            let dy = (res.ca[1] - protein.ca[[1, i]]).abs();
            let dz = (res.ca[2] - protein.ca[[2, i]]).abs();
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            assert!(
                dist < tol,
                "{}: CA[{}] mismatch dist={:.3}",
                dat_code, i, dist
            );
        }
    }
}

/// Test that PDB reader works for the 18-struct expanded corpus.
#[test]
fn test_pdb_reader_18struct_nres_and_ca() {
    let cases = [
        ("1a7s", "A", "1a7sA"),
        ("1aiw", "A", "1aiwA"),
        ("1a25", "A", "1a25A"),
        ("1a12", "A", "1a12A"),
        ("1f3u", "A", "1f3uA"),
        ("1a04", "A", "1a04A"),
        ("1bbh", "A", "1bbhA"),
        ("1b3q", "A", "1b3qA"),
        ("1a17", "A", "1a17A"),
        ("1miw", "A", "1miwA"),
        ("1aop", "A", "1aopA"),
        ("1a8l", "A", "1a8lA"),
        ("1a6q", "A", "1a6qA"),
        ("1a06", "A", "1a06A"),
        ("1b3o", "B", "1b3oB"),
        ("1a0c", "A", "1a0cA"),
        ("1a4i", "A", "1a4iA"),
        ("1bco", "A", "1bcoA"),
    ];

    let mut passed = 0;
    for (pdb_code, chain, dat_code) in &cases {
        let pdb_path = format!("{}/pdb{}.ent.gz", PDB_DIR_18, pdb_code);
        let dat_path = format!("{}/{}.dat", DAT_DIR_18, dat_code);

        let backbone = match read_pdb(&pdb_path, chain, pdb_code) {
            Ok(b) => b,
            Err(e) => {
                eprintln!("  SKIP {}: PDB read error: {}", dat_code, e);
                continue;
            }
        };
        let protein = match read_dat(&dat_path) {
            Ok(p) => p,
            Err(e) => {
                eprintln!("  SKIP {}: DAT read error: {}", dat_code, e);
                continue;
            }
        };

        let nres_match = backbone.residues.len() == protein.nres;
        let tol = 0.15;
        let mut ca_ok = nres_match;
        if nres_match {
            for (i, res) in backbone.residues.iter().enumerate() {
                let dx = (res.ca[0] - protein.ca[[0, i]]).abs();
                let dy = (res.ca[1] - protein.ca[[1, i]]).abs();
                let dz = (res.ca[2] - protein.ca[[2, i]]).abs();
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                if dist >= tol {
                    ca_ok = false;
                    break;
                }
            }
        }
        if nres_match && ca_ok {
            passed += 1;
        }
    }

    eprintln!("PDB reader 18-struct: {}/18", passed);
    assert!(passed >= 16, "PDB reader: only {}/18 passed", passed);
}

/// Full import pipeline validation for 5-struct corpus.
///
/// For each structure, imports from PDB via the full pipeline (DSSP → secstr → domain),
/// then compares against the reference .dat file. Checks:
/// - nres (must match exactly)
/// - CA coordinates (within 0.15 Å)
/// - nseg, na, nb (report match rate, don't require exact match since DSSP may differ)
/// - segment ranges (report overlap)
/// - secstr types (report match rate)
#[test]
fn test_import_pipeline_5struct() {
    let cases = [
        ("101m", "A", "101mA"),
        ("1a00", "A", "1a00A"),
        ("1a87", "A", "1a87A"),
        ("1all", "A", "1allA"),
        ("1bin", "A", "1binA"),
    ];

    eprintln!("\n=== Import Pipeline Validation (5-struct) ===");

    let mut total_seg_match = 0;
    let mut total_seg_ref = 0;

    for (pdb_code, chain, dat_code) in &cases {
        let pdb_path = format!("{}/pdb{}.ent.gz", PDB_DIR_5, pdb_code);
        let dat_path = format!("{}/{}.dat", DAT_DIR_5, dat_code);

        let imported = import_pdb(&pdb_path, chain, pdb_code).unwrap_or_else(|e| {
            panic!("Failed to import {}: {}", dat_code, e);
        });
        let reference = read_dat(&dat_path).unwrap_or_else(|e| {
            panic!("Failed to read DAT {}: {}", dat_code, e);
        });

        // nres must match exactly
        assert_eq!(
            imported.nres, reference.nres,
            "{}: nres mismatch (imported={}, ref={})",
            dat_code, imported.nres, reference.nres
        );

        // CA coordinates must match (same PDB reader, should be exact)
        let tol = 0.15;
        for i in 0..imported.nres {
            let dx = (imported.ca[[0, i]] - reference.ca[[0, i]]).abs();
            let dy = (imported.ca[[1, i]] - reference.ca[[1, i]]).abs();
            let dz = (imported.ca[[2, i]] - reference.ca[[2, i]]).abs();
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            assert!(dist < tol, "{}: CA[{}] mismatch dist={:.3}", dat_code, i, dist);
        }

        // Compare SSE segments: count how many reference segments are "covered"
        // by an imported segment of the same type (within ±2 residues)
        let seg_tol = 2i32;
        let mut seg_matches = 0;
        for ref_seg in &reference.segments {
            let ref_start = ref_seg.start as i32;
            let ref_end = ref_seg.end as i32;
            let ref_type = ref_seg.sse_type;

            let matched = imported.segments.iter().any(|imp_seg| {
                imp_seg.sse_type == ref_type
                    && (imp_seg.start as i32 - ref_start).abs() <= seg_tol
                    && (imp_seg.end as i32 - ref_end).abs() <= seg_tol
            });
            if matched {
                seg_matches += 1;
            }
        }

        let ref_nseg = reference.nseg;
        total_seg_match += seg_matches;
        total_seg_ref += ref_nseg;

        let secstr_str_imp: String = imported.secstr.iter().map(|s| s.to_char()).collect();
        let secstr_str_ref: String = reference.secstr.iter().map(|s| s.to_char()).collect();

        eprintln!(
            "  {}: nres={} nseg imp/ref={}/{} na={}/{} nb={}/{} seg_match={}/{} secstr imp={} ref={}",
            dat_code,
            imported.nres,
            imported.nseg, reference.nseg,
            imported.na, reference.na,
            imported.nb, reference.nb,
            seg_matches, ref_nseg,
            secstr_str_imp, secstr_str_ref,
        );

        // Domain tree should have at least 1 node
        assert!(
            !imported.domain_tree.is_empty(),
            "{}: domain tree is empty",
            dat_code
        );
    }

    let seg_pct = if total_seg_ref > 0 {
        100.0 * total_seg_match as f64 / total_seg_ref as f64
    } else {
        0.0
    };
    eprintln!(
        "\nSegment overlap: {}/{} ({:.1}%)",
        total_seg_match, total_seg_ref, seg_pct
    );

    // DSSP reimplementation may not match DaliLite's dsspcmbi exactly,
    // but should capture most SSEs. Require at least 60% segment overlap.
    assert!(
        seg_pct >= 60.0,
        "Segment overlap too low: {:.1}% (need >= 60%)",
        seg_pct
    );
}

/// Full import pipeline validation for 18-struct expanded corpus.
#[test]
fn test_import_pipeline_18struct() {
    let cases = [
        ("1a7s", "A", "1a7sA"),
        ("1aiw", "A", "1aiwA"),
        ("1a25", "A", "1a25A"),
        ("1a12", "A", "1a12A"),
        ("1f3u", "A", "1f3uA"),
        ("1a04", "A", "1a04A"),
        ("1bbh", "A", "1bbhA"),
        ("1b3q", "A", "1b3qA"),
        ("1a17", "A", "1a17A"),
        ("1miw", "A", "1miwA"),
        ("1aop", "A", "1aopA"),
        ("1a8l", "A", "1a8lA"),
        ("1a6q", "A", "1a6qA"),
        ("1a06", "A", "1a06A"),
        ("1b3o", "B", "1b3oB"),
        ("1a0c", "A", "1a0cA"),
        ("1a4i", "A", "1a4iA"),
        ("1bco", "A", "1bcoA"),
    ];

    eprintln!("\n=== Import Pipeline Validation (18-struct) ===");

    let mut structures_ok = 0;
    let mut total_seg_match = 0;
    let mut total_seg_ref = 0;

    for (pdb_code, chain, dat_code) in &cases {
        let pdb_path = format!("{}/pdb{}.ent.gz", PDB_DIR_18, pdb_code);
        let dat_path = format!("{}/{}.dat", DAT_DIR_18, dat_code);

        let imported = match import_pdb(&pdb_path, chain, pdb_code) {
            Ok(p) => p,
            Err(e) => {
                eprintln!("  SKIP {}: import error: {}", dat_code, e);
                continue;
            }
        };
        let reference = match read_dat(&dat_path) {
            Ok(p) => p,
            Err(e) => {
                eprintln!("  SKIP {}: DAT read error: {}", dat_code, e);
                continue;
            }
        };

        let nres_ok = imported.nres == reference.nres;

        // Segment overlap
        let seg_tol = 2i32;
        let mut seg_matches = 0;
        if nres_ok {
            for ref_seg in &reference.segments {
                let ref_start = ref_seg.start as i32;
                let ref_end = ref_seg.end as i32;
                let ref_type = ref_seg.sse_type;

                let matched = imported.segments.iter().any(|imp_seg| {
                    imp_seg.sse_type == ref_type
                        && (imp_seg.start as i32 - ref_start).abs() <= seg_tol
                        && (imp_seg.end as i32 - ref_end).abs() <= seg_tol
                });
                if matched {
                    seg_matches += 1;
                }
            }
            total_seg_match += seg_matches;
            total_seg_ref += reference.nseg;
        }

        let status = if nres_ok { "OK" } else { "NRES_MISMATCH" };
        if nres_ok {
            structures_ok += 1;
        }

        eprintln!(
            "  {}: {} nres={}/{} nseg={}/{} seg_match={}/{}",
            dat_code, status,
            imported.nres, reference.nres,
            imported.nseg, reference.nseg,
            seg_matches, reference.nseg,
        );
    }

    let seg_pct = if total_seg_ref > 0 {
        100.0 * total_seg_match as f64 / total_seg_ref as f64
    } else {
        0.0
    };
    eprintln!(
        "\nStructures with matching nres: {}/18",
        structures_ok
    );
    eprintln!(
        "Segment overlap: {}/{} ({:.1}%)",
        total_seg_match, total_seg_ref, seg_pct
    );

    assert!(structures_ok >= 16, "Only {}/18 structures matched nres", structures_ok);
    assert!(seg_pct >= 50.0, "Segment overlap too low: {:.1}%", seg_pct);
}

/// Test write_dat round-trip: read .dat → write → read back → compare.
#[test]
fn test_write_dat_roundtrip() {
    let cases = ["101mA", "1a00A", "1a87A", "1allA", "1binA"];

    for code in &cases {
        let dat_path = format!("{}/{}.dat", DAT_DIR_5, code);
        let original = read_dat(&dat_path).unwrap_or_else(|e| {
            panic!("Failed to read {}: {}", code, e);
        });

        // Write to temp file
        let tmp_path = format!("/tmp/test_write_dat_{}.dat", code);
        write_dat(&original, &tmp_path).unwrap_or_else(|e| {
            panic!("Failed to write {}: {}", code, e);
        });

        // Read back
        let roundtrip = read_dat(&tmp_path).unwrap_or_else(|e| {
            panic!("Failed to read back {}: {}", code, e);
        });

        // Verify fields match
        assert_eq!(roundtrip.code, original.code, "{}: code mismatch", code);
        assert_eq!(roundtrip.nres, original.nres, "{}: nres mismatch", code);
        assert_eq!(roundtrip.nseg, original.nseg, "{}: nseg mismatch", code);
        assert_eq!(roundtrip.na, original.na, "{}: na mismatch", code);
        assert_eq!(roundtrip.nb, original.nb, "{}: nb mismatch", code);

        // CA coordinates within tolerance (f8.1 format loses precision)
        for j in 0..original.nres {
            for d in 0..3 {
                let diff = (roundtrip.ca[[d, j]] - original.ca[[d, j]]).abs();
                assert!(diff < 0.1, "{}: CA[{},{}] diff={:.3}", code, d, j, diff);
            }
        }

        // Domain tree
        assert_eq!(roundtrip.domain_tree.len(), original.domain_tree.len(),
                   "{}: domain tree length mismatch", code);
        for (i, (rt, orig)) in roundtrip.domain_tree.iter()
            .zip(original.domain_tree.iter()).enumerate()
        {
            assert_eq!(rt.nseg, orig.nseg,
                       "{}: domain[{}] nseg mismatch ({} vs {})", code, i, rt.nseg, orig.nseg);
            assert_eq!(rt.segments, orig.segments,
                       "{}: domain[{}] segments mismatch", code, i);
        }

        // Clean up
        let _ = std::fs::remove_file(&tmp_path);
        eprintln!("  {}: write_dat roundtrip OK", code);
    }
}
