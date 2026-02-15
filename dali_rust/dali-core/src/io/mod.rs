pub mod dat;
pub mod domain;
pub mod dssp;
pub mod import;
pub mod pdb;
pub mod secstr;

pub use dat::{read_dat, write_dat, DatError};
pub use domain::build_domain_tree;
pub use dssp::compute_dssp;
pub use import::{import_pdb, ImportError};
pub use pdb::{read_pdb, BackboneData, BackboneResidue, PdbReadError};
pub use secstr::assign_segments;
