//! DSSP secondary structure assignment (Kabsch-Sander algorithm).
//!
//! Reference: Kabsch W, Sander C. "Dictionary of protein secondary structure:
//! pattern recognition of hydrogen-bonded and geometrical features."
//! Biopolymers 22:2577-2637, 1983.
//!
//! Assigns 8-state secondary structure per residue: H, G, I, E, B, T, S, or ' '.

use crate::io::pdb::BackboneResidue;

/// Coulomb constant Q = q1*q2*f = 0.42 * 0.20 * 332 * 1000 cal/mol = 27888 cal/mol
/// (positive: the K-S formula E = Q*(1/rON + 1/rCH - 1/rOH - 1/rCN) gives E < 0 for H-bonds)
const Q: f64 = 27888.0;
/// H-bond energy cutoff (cal/mol). H-bond if E < HBHIGH.
const HBHIGH: f64 = -500.0;
/// Minimum distance for H-bond calculation (Å).
const MIN_DIST: f64 = 0.5;
/// Maximum CA-CA distance for H-bond search (Å).
const CA_DIST_CUTOFF: f64 = 9.0;

/// An H-bond record: (partner_index, energy).
/// partner_index is the residue index (0-based), energy in cal/mol.
#[derive(Debug, Clone, Copy)]
struct HBond {
    partner: i32, // -1 = no bond
    energy: f64,
}

impl HBond {
    fn none() -> Self {
        HBond {
            partner: -1,
            energy: 0.0,
        }
    }
}

/// Distance between two 3D points.
fn dist(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Estimate NH hydrogen position from backbone geometry.
///
/// Standard DSSP formula: H = N + unit(prev_C - prev_O)
/// Places H 1.0 Å from N in the direction opposite to the C=O bond
/// of the preceding residue (trans peptide bond geometry).
fn estimate_h_position(n: &[f64; 3], prev_c: &[f64; 3], prev_o: &[f64; 3]) -> [f64; 3] {
    // H = N + normalize(prev_C - prev_O)
    // This is the standard mkdssp formula.
    let co_x = prev_c[0] - prev_o[0];
    let co_y = prev_c[1] - prev_o[1];
    let co_z = prev_c[2] - prev_o[2];
    let co_len = (co_x * co_x + co_y * co_y + co_z * co_z).sqrt();
    if co_len < 0.001 {
        return *n; // degenerate case
    }
    [
        n[0] + co_x / co_len,
        n[1] + co_y / co_len,
        n[2] + co_z / co_len,
    ]
}

/// Calculate H-bond energy between donor (NH at residue i) and acceptor (CO at residue j).
///
/// E = Q * (1/r_ON + 1/r_CH - 1/r_OH - 1/r_CN)
///
/// Returns energy in cal/mol, or None if atoms are too close.
fn hbond_energy(
    h: &[f64; 3],  // H position (donor)
    n: &[f64; 3],  // N position (donor)
    c: &[f64; 3],  // C position (acceptor)
    o: &[f64; 3],  // O position (acceptor)
) -> Option<f64> {
    let r_on = dist(o, n);
    let r_ch = dist(c, h);
    let r_oh = dist(o, h);
    let r_cn = dist(c, n);

    if r_on < MIN_DIST || r_ch < MIN_DIST || r_oh < MIN_DIST || r_cn < MIN_DIST {
        return None;
    }

    let energy = Q * (1.0 / r_on + 1.0 / r_ch - 1.0 / r_oh - 1.0 / r_cn);
    Some(energy)
}

/// Bridge type between two residues.
#[derive(Debug, Clone, Copy, PartialEq)]
enum BridgeType {
    None,
    Parallel,
    AntiParallel,
}

/// Test if there's an H-bond from donor i to acceptor j.
fn test_bond(acc: &[[HBond; 2]], i: usize, j: usize) -> bool {
    // Check if residue i donates to j (i's NH → j's CO)
    let nres = acc.len();
    if i >= nres || j >= nres {
        return false;
    }
    acc[i][0].partner == j as i32 || acc[i][1].partner == j as i32
}

/// Test for a bridge between residues i and j.
///
/// Uses the donor array only: donor[i] = best 2 acceptors for NH_i.
/// Hbond(a,b) in K-S notation means "CO of a accepts NH of b" = "b donates to a".
fn test_bridge(
    donor: &[[HBond; 2]],
    i: usize,
    j: usize,
) -> BridgeType {
    let nres = donor.len();
    if i < 1 || j < 1 || i >= nres - 1 || j >= nres - 1 {
        return BridgeType::None;
    }

    // Parallel bridge: (HBond(i-1,j) AND HBond(j,i+1)) OR (HBond(j-1,i) AND HBond(i,j+1))
    let parallel = (test_bond(donor, i + 1, j) && test_bond(donor, j, i - 1))
        || (test_bond(donor, j + 1, i) && test_bond(donor, i, j - 1));

    // Antiparallel bridge: (HBond(i,j) AND HBond(j,i)) OR (HBond(i-1,j+1) AND HBond(j-1,i+1))
    let antiparallel = (test_bond(donor, i, j) && test_bond(donor, j, i))
        || (test_bond(donor, i + 1, j - 1) && test_bond(donor, j + 1, i - 1));

    if parallel {
        BridgeType::Parallel
    } else if antiparallel {
        BridgeType::AntiParallel
    } else {
        BridgeType::None
    }
}

/// Compute DSSP 8-state secondary structure from backbone coordinates.
///
/// Input: slice of BackboneResidue with N, CA, C, O coordinates.
/// Output: Vec<u8> of DSSP codes, one per residue.
pub fn compute_dssp(residues: &[BackboneResidue]) -> Vec<u8> {
    let nres = residues.len();
    if nres < 2 {
        return vec![b' '; nres];
    }

    // Step 1: Estimate H positions
    let mut h_pos: Vec<[f64; 3]> = Vec::with_capacity(nres);
    h_pos.push(residues[0].n); // first residue has no preceding C=O
    for i in 1..nres {
        // Check for chain break: if CA-CA distance > 4.5 Å, no H-bond
        let ca_dist = dist(&residues[i].ca, &residues[i - 1].ca);
        if ca_dist > 5.5 {
            h_pos.push(residues[i].n); // can't estimate H across chain break
        } else {
            h_pos.push(estimate_h_position(
                &residues[i].n,
                &residues[i - 1].c,
                &residues[i - 1].o,
            ));
        }
    }

    // Step 2: Calculate H-bonds
    // donor[i] = best 2 acceptors for NH_i (residue i is the donor)
    let mut donor: Vec<[HBond; 2]> = vec![[HBond::none(); 2]; nres];

    for i in 1..nres {
        // Proline cannot donate (no NH)
        if residues[i].resname == "PRO" {
            continue;
        }

        for j in 0..nres {
            if i == j {
                continue;
            }
            // Skip if too far apart
            let ca_d = dist(&residues[i].ca, &residues[j].ca);
            if ca_d > CA_DIST_CUTOFF {
                continue;
            }
            // Skip if |i-j| < 3 (too close in sequence)
            if (i as i64 - j as i64).unsigned_abs() < 3 {
                continue;
            }

            // Calculate energy: i donates NH to j's CO
            if let Some(energy) = hbond_energy(
                &h_pos[i],
                &residues[i].n,
                &residues[j].c,
                &residues[j].o,
            ) {
                if energy < HBHIGH {
                    // Insert into best-2 list for donor i
                    if energy < donor[i][0].energy || donor[i][0].partner == -1 {
                        donor[i][1] = donor[i][0];
                        donor[i][0] = HBond {
                            partner: j as i32,
                            energy,
                        };
                    } else if energy < donor[i][1].energy || donor[i][1].partner == -1 {
                        donor[i][1] = HBond {
                            partner: j as i32,
                            energy,
                        };
                    }
                }
            }
        }
    }

    // DSSP priority: H > B > E > G > I > T > S > ' '
    // Assignment order: bridges/sheets first, then helices, then turns, then bends.

    let mut ss = vec![b' '; nres];

    // Step 3: Detect bridges and sheets (highest non-helix priority)
    let mut bridge_partner1: Vec<i32> = vec![-1; nres];
    let mut bridge_partner2: Vec<i32> = vec![-1; nres];

    for i in 1..nres - 1 {
        for j in (i + 2)..nres - 1 {
            let bt = test_bridge(&donor, i, j);
            if bt != BridgeType::None {
                if bridge_partner1[i] == -1 {
                    bridge_partner1[i] = j as i32;
                } else if bridge_partner2[i] == -1 {
                    bridge_partner2[i] = j as i32;
                }
                if bridge_partner1[j] == -1 {
                    bridge_partner1[j] = i as i32;
                } else if bridge_partner2[j] == -1 {
                    bridge_partner2[j] = i as i32;
                }
            }
        }
    }

    // Assign E (extended strand) and B (isolated bridge)
    for i in 1..nres - 1 {
        if bridge_partner1[i] == -1 {
            continue;
        }
        let is_ladder = {
            let prev_has = i > 0 && bridge_partner1[i - 1] != -1;
            let next_has = i < nres - 1 && bridge_partner1[i + 1] != -1;
            prev_has || next_has
        };
        if is_ladder {
            ss[i] = b'E';
        } else {
            ss[i] = b'B';
        }
    }

    // Step 4: Detect turns and helices
    // n-turn detection: residue i has n-turn if i+n donates to i
    let mut has_turn = vec![[false; 3]; nres]; // [3-turn, 4-turn, 5-turn]

    for n_idx in 0..3usize {
        let n = n_idx + 3; // 3, 4, 5
        for i in 1..nres.saturating_sub(n) {
            if test_bond(&donor, i + n, i) {
                has_turn[i][n_idx] = true;
            }
        }
    }

    // Helix assignment: consecutive n-turns → helix
    // Process in priority order: H (4-turn) first, then G (3-turn), then I (5-turn)
    // H overrides everything; G overrides I,T,S,' '; I overrides T,S,' '
    let helix_order = [(1usize, b'H'), (0usize, b'G'), (2usize, b'I')];

    for &(n_idx, code) in &helix_order {
        let n = n_idx + 3;
        for i in 1..nres.saturating_sub(n) {
            if i > 0 && has_turn[i][n_idx] && has_turn[i - 1][n_idx] {
                for k in i..=(i + n - 1).min(nres - 1) {
                    let cur = ss[k];
                    let dominated = match code {
                        b'H' => true, // H overrides everything
                        b'G' => cur == b'I' || cur == b'T' || cur == b'S' || cur == b' ',
                        b'I' => cur == b'T' || cur == b'S' || cur == b' ',
                        _ => false,
                    };
                    if dominated {
                        ss[k] = code;
                    }
                }
            }
        }
    }

    // Mark turns (T): only if current is S or ' '
    for n_idx in 0..3usize {
        let n = n_idx + 3;
        for i in 1..nres.saturating_sub(n) {
            if has_turn[i][n_idx] {
                for k in (i + 1)..=(i + n - 1).min(nres - 1) {
                    if ss[k] == b' ' || ss[k] == b'S' {
                        ss[k] = b'T';
                    }
                }
            }
        }
    }

    // Step 5: Detect bends (S)
    // A bend occurs at residue i if the angle between CA(i-2)→CA(i) and CA(i)→CA(i+2) > 70°
    for i in 2..nres.saturating_sub(2) {
        if ss[i] != b' ' {
            continue;
        }
        let v1 = [
            residues[i].ca[0] - residues[i - 2].ca[0],
            residues[i].ca[1] - residues[i - 2].ca[1],
            residues[i].ca[2] - residues[i - 2].ca[2],
        ];
        let v2 = [
            residues[i + 2].ca[0] - residues[i].ca[0],
            residues[i + 2].ca[1] - residues[i].ca[1],
            residues[i + 2].ca[2] - residues[i].ca[2],
        ];
        let dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
        let len1 = (v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]).sqrt();
        let len2 = (v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]).sqrt();
        if len1 > 0.001 && len2 > 0.001 {
            let cos_angle = dot / (len1 * len2);
            let angle = cos_angle.acos().to_degrees();
            if angle > 70.0 {
                ss[i] = b'S';
            }
        }
    }

    ss
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hbond_energy_basic() {
        // Place atoms in a configuration that should give a weak H-bond
        let h = [0.0, 0.0, 0.0];
        let n = [-1.0, 0.0, 0.0];
        let c = [3.0, 0.0, 0.0];
        let o = [4.0, 0.0, 0.0];
        let e = hbond_energy(&h, &n, &c, &o);
        assert!(e.is_some());
    }

    #[test]
    fn test_dssp_empty() {
        let result = compute_dssp(&[]);
        assert!(result.is_empty());
    }

    #[test]
    fn test_dssp_on_real_pdb() {
        let pdb_path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../..",
            "/validation/fixtures/structures/pdb/pdb101m.ent.gz"
        );
        let backbone = crate::io::pdb::read_pdb(pdb_path, "A", "101m").unwrap();
        eprintln!("Loaded {} residues", backbone.residues.len());

        // Check backbone geometry
        for i in 0..3.min(backbone.residues.len()) {
            let r = &backbone.residues[i];
            eprintln!("  res[{}] {}: N=({:.1},{:.1},{:.1}) CA=({:.1},{:.1},{:.1}) C=({:.1},{:.1},{:.1}) O=({:.1},{:.1},{:.1})",
                i, r.resname,
                r.n[0], r.n[1], r.n[2],
                r.ca[0], r.ca[1], r.ca[2],
                r.c[0], r.c[1], r.c[2],
                r.o[0], r.o[1], r.o[2]);
        }

        let dssp = compute_dssp(&backbone.residues);
        let dssp_str: String = dssp.iter().map(|&b| b as char).collect();
        eprintln!("DSSP (len={}): '{}'", dssp.len(), dssp_str);

        // Count types
        let h_count = dssp.iter().filter(|&&b| b == b'H').count();
        let e_count = dssp.iter().filter(|&&b| b == b'E').count();
        let g_count = dssp.iter().filter(|&&b| b == b'G').count();
        let t_count = dssp.iter().filter(|&&b| b == b'T').count();
        let s_count = dssp.iter().filter(|&&b| b == b'S').count();
        let space_count = dssp.iter().filter(|&&b| b == b' ').count();
        eprintln!("  H={} E={} G={} T={} S={} ' '={}", h_count, e_count, g_count, t_count, s_count, space_count);

        // 101m is myoglobin, an all-alpha protein. Should have significant helix content.
        assert!(h_count > 50, "Expected >50 H residues for myoglobin, got {}", h_count);
    }

    #[test]
    fn test_dssp_beta_sheets() {
        // 1a87A is a TIM barrel (α/β protein) — should detect both helices and strands
        let pdb_path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../..",
            "/validation/fixtures/structures/pdb/pdb1a87.ent.gz"
        );
        let backbone = crate::io::pdb::read_pdb(pdb_path, "A", "1a87").unwrap();
        let dssp = compute_dssp(&backbone.residues);

        let h_count = dssp.iter().filter(|&&b| b == b'H').count();
        let e_count = dssp.iter().filter(|&&b| b == b'E').count();

        // TIM barrel should have both alpha helices and beta strands
        assert!(h_count > 30, "Expected >30 H residues for TIM barrel, got {}", h_count);
        assert!(e_count > 10, "Expected >10 E residues for TIM barrel, got {}", e_count);
    }
}
