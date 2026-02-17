pub mod protein;
pub mod alignment;
pub mod distance;

pub use protein::{DomainNode, NodeType, Protein, ResidNumbering, Segment, SseType};
pub use alignment::{
    Alignment, AlignmentBlock, DccpEntry, ParsiHit, ScoredAlignment, compress_blocks,
};
pub use distance::{DaliconDistMat, DistMatScale10, DistMatScale100, ParsiDistMat};

use ndarray::{Array1, Array2};

/// Result of a single hit from iterative database search.
#[derive(Debug, Clone)]
pub struct SearchHit {
    pub cd2: String,
    pub zscore: f64,
    pub score: f64,
    pub rmsd: f64,
    pub nblock: usize,
    pub blocks: Vec<AlignmentBlock>,
    pub rotation: Array2<f64>,
    pub translation: Array1<f64>,
    pub alignments: Vec<(usize, usize)>,  // (query_1based, template_1based)
    pub round: usize,
}

/// Minimum helix length for segment reconstruction (matches io/secstr.rs LENH).
const MASK_LENH: usize = 8;
/// Minimum strand length for segment reconstruction (matches io/secstr.rs LENE).
const MASK_LENE: usize = 6;

/// Create a new Protein containing only the specified residues.
///
/// `keep_indices` are 0-based indices into the original protein's residue array.
/// Reconstructs segments from the per-residue secondary structure, applies
/// minimum-length filters, and creates a single-node domain tree.
pub fn mask_protein(
    protein: &Protein,
    keep_indices: &[usize],
    new_code: &str,
) -> Protein {
    let new_nres = keep_indices.len();

    // Subset CA coordinates
    let mut ca = Array2::<f64>::zeros((3, new_nres));
    for (new_i, &old_i) in keep_indices.iter().enumerate() {
        for d in 0..3 {
            ca[[d, new_i]] = protein.ca[[d, old_i]];
        }
    }

    // Reconstruct per-residue secstr from original segments
    let mut residue_ss = vec![b'L'; protein.nres];
    for seg in &protein.segments {
        let c = seg.sse_type.to_char() as u8;
        for i in seg.start..=seg.end {
            residue_ss[i as usize - 1] = c; // segments are 1-based
        }
    }

    // Subset to kept residues
    let subset_ss: Vec<u8> = keep_indices.iter().map(|&i| residue_ss[i]).collect();

    // Scan for H/E runs, filter by min length → new Segment list
    let mut segments = Vec::new();
    let mut secstr = Vec::new();
    if !subset_ss.is_empty() {
        let mut run_start = 0;
        let mut run_char = subset_ss[0];
        for i in 1..new_nres {
            if subset_ss[i] != run_char {
                if run_char == b'H' || run_char == b'E' {
                    let len = i - run_start;
                    let min_len = if run_char == b'H' { MASK_LENH } else { MASK_LENE };
                    if len >= min_len {
                        let sse_type = if run_char == b'H' { SseType::Helix } else { SseType::Strand };
                        secstr.push(sse_type);
                        segments.push(Segment {
                            start: (run_start + 1) as u32,
                            end: i as u32,
                            sse_type,
                            check_start: 0,
                            check_end: 0,
                            checkx: 0,
                        });
                    }
                }
                run_start = i;
                run_char = subset_ss[i];
            }
        }
        // Final run
        if run_char == b'H' || run_char == b'E' {
            let len = new_nres - run_start;
            let min_len = if run_char == b'H' { MASK_LENH } else { MASK_LENE };
            if len >= min_len {
                let sse_type = if run_char == b'H' { SseType::Helix } else { SseType::Strand };
                secstr.push(sse_type);
                segments.push(Segment {
                    start: (run_start + 1) as u32,
                    end: new_nres as u32,
                    sse_type,
                    check_start: 0,
                    check_end: 0,
                    checkx: 0,
                });
            }
        }
    }

    // Compute check ranges (same logic as io/secstr.rs)
    for (i, seg) in segments.iter_mut().enumerate() {
        let len = (seg.end - seg.start + 1) as usize;
        let x = if secstr[i] == SseType::Helix {
            if len >= 9 { 0 } else if len >= 7 { 1 } else { 2 }
        } else {
            if len >= 5 { 0 } else if len >= 3 { 1 } else { 2 }
        };
        seg.checkx = x;
        seg.check_start = (seg.start as i32 - x).max(1) as u32;
        seg.check_end = (seg.end as i32 + x).min(new_nres as i32) as u32;
    }

    let na = secstr.iter().filter(|&&t| t == SseType::Helix).count();
    let nb = secstr.iter().filter(|&&t| t == SseType::Strand).count();
    let nseg = segments.len();

    // Subset resid_map
    let resid_map: Vec<i32> = keep_indices.iter().map(|&i| protein.resid_map[i]).collect();

    // Subset sequence
    let sequence = if protein.sequence.is_empty() {
        String::new()
    } else {
        let chars: Vec<char> = protein.sequence.chars().collect();
        keep_indices.iter().map(|&i| chars[i]).collect()
    };

    // Single-node domain tree covering all residues
    let domain_tree = vec![DomainNode {
        index: 1,
        node_type: NodeType::Root,
        left_child: 0,
        right_child: 0,
        nseg: 1,
        segments: vec![(1, new_nres as u32)],
    }];

    Protein {
        code: new_code.to_string(),
        nres: new_nres,
        nseg,
        na,
        nb,
        segments,
        secstr,
        ca,
        sequence,
        domain_tree,
        resid_map,
        numbering: protein.numbering,
    }
}
