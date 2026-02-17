pub mod types;
pub mod numerics;
pub mod io;
pub mod store;

pub mod wolf;
pub mod dp;
pub mod dalicon;
pub mod parsi;
pub mod filter95;
pub mod pipeline;

// Re-export key types
pub use types::{
    Protein, Segment, SseType, DomainNode, NodeType, ResidNumbering,
    Alignment, AlignmentBlock, ScoredAlignment, DccpEntry, ParsiHit,
    DistMatScale10, DistMatScale100, DaliconDistMat, ParsiDistMat,
    compress_blocks, SearchHit, mask_protein,
};
pub use store::ProteinStore;
