/// A single alignment block mapping a range in protein 1 to a range in protein 2.
#[derive(Debug, Clone, PartialEq)]
pub struct AlignmentBlock {
    pub l1: u32, // start in protein 1 (1-based)
    pub r1: u32, // end in protein 1 (1-based)
    pub l2: u32, // start in protein 2 (1-based)
    pub r2: u32, // end in protein 2 (1-based)
}

/// An alignment between two proteins (no score information).
#[derive(Debug, Clone)]
pub struct Alignment {
    pub cd1: String,
    pub cd2: String,
    pub blocks: Vec<AlignmentBlock>,
}

/// An alignment with Z-score and RMSD from DP scoring.
#[derive(Debug, Clone)]
pub struct ScoredAlignment {
    pub alignment: Alignment,
    pub score: f64,
    pub zscore: f64,
    pub rmsd: f64,
}

/// A hit from PARSI exhaustive search.
#[derive(Debug, Clone)]
pub struct ParsiHit {
    pub cd1: String,
    pub cd2: String,
    pub idom: usize,
    pub score: i32,
    pub ranges_cd1: Vec<(i32, i32)>, // (start, end) per segment
    pub ranges_cd2: Vec<(i32, i32)>,
}

/// Final DCCP output entry.
#[derive(Debug, Clone)]
pub struct DccpEntry {
    pub cd1: String,
    pub cd2: String,
    pub score: f64,
    pub zscore: f64,
    pub rmsd: f64,
    pub blocks: Vec<AlignmentBlock>,
}

/// Merge consecutive single-residue alignment blocks.
///
/// Exact translation of comparemodules.f compressblocks().
pub fn compress_blocks(blocks: &[AlignmentBlock]) -> Vec<AlignmentBlock> {
    if blocks.is_empty() {
        return Vec::new();
    }

    let nblock = blocks.len();
    let mut result: Vec<AlignmentBlock> = Vec::new();
    let mut i = 0;

    while i < nblock {
        let mut block = blocks[i].clone();

        if i < nblock - 1 {
            // Try to merge consecutive single-residue blocks
            while blocks[i].l1 == blocks[i].r1
                && blocks[i].l2 == blocks[i].r2
                && blocks[i + 1].l1 == blocks[i].r1 + 1
                && blocks[i + 1].l2 == blocks[i].r2 + 1
            {
                if i + 1 >= nblock - 1 {
                    break;
                }
                i += 1;
            }
            block.r1 = blocks[i].r1;
            block.r2 = blocks[i].r2;
        }

        result.push(block);
        i += 1;
    }

    result
}
