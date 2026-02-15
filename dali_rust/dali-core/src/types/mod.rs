pub mod protein;
pub mod alignment;
pub mod distance;

pub use protein::{DomainNode, NodeType, Protein, Segment, SseType};
pub use alignment::{
    Alignment, AlignmentBlock, DccpEntry, ParsiHit, ScoredAlignment, compress_blocks,
};
pub use distance::{DaliconDistMat, DistMatScale10, DistMatScale100, ParsiDistMat};
