use ndarray::Array2;

/// Secondary structure element type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SseType {
    Helix,
    Strand,
}

impl SseType {
    /// Parse from single-character DSSP code used in .dat files.
    pub fn from_char(c: char) -> Option<Self> {
        match c {
            'H' => Some(SseType::Helix),
            'E' => Some(SseType::Strand),
            _ => None,
        }
    }

    pub fn to_char(self) -> char {
        match self {
            SseType::Helix => 'H',
            SseType::Strand => 'E',
        }
    }
}

/// A secondary structure segment (SSE).
#[derive(Debug, Clone)]
pub struct Segment {
    pub start: u32,       // 1-based inclusive
    pub end: u32,         // 1-based inclusive
    pub sse_type: SseType,
    pub check_start: u32, // from .dat file
    pub check_end: u32,
    pub checkx: i32,
}

/// Node type in the domain decomposition tree.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NodeType {
    Root,  // '*'
    Split, // '+'
    Minus, // '-'
}

impl NodeType {
    pub fn from_char(c: char) -> Self {
        match c {
            '*' => NodeType::Root,
            '+' => NodeType::Split,
            _ => NodeType::Minus,
        }
    }

    pub fn to_char(self) -> char {
        match self {
            NodeType::Root => '*',
            NodeType::Split => '+',
            NodeType::Minus => '-',
        }
    }
}

/// A node in the hierarchical domain decomposition tree.
#[derive(Debug, Clone)]
pub struct DomainNode {
    pub index: usize,            // 1-based node index
    pub node_type: NodeType,
    pub left_child: usize,       // 0 = leaf
    pub right_child: usize,      // 0 = leaf
    pub nseg: usize,             // number of segments in this domain
    pub segments: Vec<(u32, u32)>, // (start, end) pairs, 1-based
}

/// How residues in `resid_map` are numbered.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ResidNumbering {
    /// Sequential 1..=nres (from .dat files, which don't store PDB numbering).
    Sequential,
    /// PDB residue serial numbers (from `import_pdb`, may have gaps or not start at 1).
    Pdb,
}

/// A protein structure with all data needed by the DALI pipeline.
#[derive(Debug, Clone)]
pub struct Protein {
    pub code: String,                // e.g. "101mA"
    pub nres: usize,
    pub nseg: usize,                 // number of SSEs
    pub na: usize,                   // number of helices
    pub nb: usize,                   // number of strands
    pub segments: Vec<Segment>,      // nseg SSEs
    pub secstr: Vec<SseType>,        // one per SSE
    pub ca: Array2<f64>,             // (3, nres) CA coordinates
    pub sequence: String,            // amino acid sequence (if available)
    pub domain_tree: Vec<DomainNode>, // hierarchical decomposition
    pub resid_map: Vec<i32>,         // residue identifiers (see numbering field)
    pub numbering: ResidNumbering,   // how resid_map values should be interpreted
}

impl Protein {
    /// Get filtered domain definitions for DP scoring.
    ///
    /// Returns domains where node index == 1 or node_type is Root/Split.
    /// This matches the Fortran filter: `if(i.eq.1.or.node_type(i).eq.'+'.or.node_type(i).eq.'*')`
    pub fn filtered_domains(&self) -> Vec<&DomainNode> {
        self.domain_tree
            .iter()
            .filter(|n| {
                n.index == 1
                    || n.node_type == NodeType::Root
                    || n.node_type == NodeType::Split
            })
            .collect()
    }
}
