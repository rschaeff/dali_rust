//! Domain decomposition from CA coordinates.
//!
//! Simplified version of the PUU (Protein Unfolding Units) algorithm.
//! Produces a hierarchical binary tree of domain decomposition based on
//! CA-CA contact density, compatible with DaliLite .dat format.
//!
//! Reference: Holm L, Sander C. "Parser for Protein Folding Units."
//! Proteins 19:256-268, 1994.

use crate::types::{DomainNode, NodeType, Segment};

/// Contact distance cutoff for CA-CA contacts (Angstroms).
const CONTACT_DIST: f64 = 8.0;
/// Minimum sequence separation for contacts.
const MIN_SEQ_SEP: usize = 8;
/// Minimum domain size (residues) to attempt splitting.
const MIN_DOMAIN_SIZE: usize = 40;

/// A contact pair: (residue_i, residue_j, weight).
struct Contact {
    i: usize,
    j: usize,
}

/// Build a contact list from CA coordinates.
///
/// Contact: CA_i - CA_j distance < CONTACT_DIST and |i-j| > MIN_SEQ_SEP.
fn build_contacts(ca: &[[f64; 3]]) -> Vec<Contact> {
    let nres = ca.len();
    let mut contacts = Vec::new();

    for i in 0..nres {
        for j in (i + MIN_SEQ_SEP)..nres {
            let dx = ca[i][0] - ca[j][0];
            let dy = ca[i][1] - ca[j][1];
            let dz = ca[i][2] - ca[j][2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            if dist < CONTACT_DIST {
                contacts.push(Contact { i, j });
            }
        }
    }
    contacts
}

/// Count contacts crossing a binary partition.
///
/// Given a set of residues split into two groups at `cut_residues`,
/// count contacts between the two groups.
///
/// Partition: residues in `set_mask` with value 1 vs value 2.
fn count_cross_contacts(contacts: &[Contact], set: &[u8]) -> usize {
    contacts
        .iter()
        .filter(|c| set[c.i] != 0 && set[c.j] != 0 && set[c.i] != set[c.j])
        .count()
}

/// Count intra-domain contacts for a set value.
fn count_intra_contacts(contacts: &[Contact], set: &[u8], val: u8) -> usize {
    contacts
        .iter()
        .filter(|c| set[c.i] == val && set[c.j] == val)
        .count()
}

/// Compute the tau score for a partition.
/// Tau = 1 - (cross_contacts / total_contacts_in_domain).
/// Higher tau = better separation.
fn compute_tau(contacts: &[Contact], set: &[u8]) -> f64 {
    let cross = count_cross_contacts(contacts, set) as f64;
    let intra1 = count_intra_contacts(contacts, set, 1) as f64;
    let intra2 = count_intra_contacts(contacts, set, 2) as f64;
    let total = cross + intra1 + intra2;

    if total < 1.0 {
        return 0.0;
    }

    // Normalized cut score: lower cross-contacts = higher tau
    let tau = (intra1 * intra2) / (total * total);
    if cross > 0.5 * total {
        0.0 // too many cross-contacts, bad split
    } else {
        tau
    }
}

/// A node in the tree being built.
struct TreeNode {
    residues: Vec<usize>,      // 0-based residue indices
    left: usize,               // child index (0 = none)
    right: usize,              // child index (0 = none)
    node_type: char,
}

/// Find the best binary cut for a set of residues.
///
/// Tries all possible sequential cut points and returns the partition
/// that maximizes the tau score.
///
/// Returns (best_tau, partition) where partition[i] = 1 or 2.
fn find_best_cut(
    residues: &[usize],
    contacts: &[Contact],
    nres: usize,
) -> (f64, Vec<u8>) {
    let n = residues.len();
    if n < 2 * MIN_DOMAIN_SIZE {
        return (0.0, vec![]);
    }

    let mut set = vec![0u8; nres]; // 0 = not in domain

    let mut best_tau = 0.0;
    let mut best_partition = vec![0u8; nres];

    // Sort residues by index
    let mut sorted = residues.to_vec();
    sorted.sort();

    // Try each sequential cut point
    for cut_pos in MIN_DOMAIN_SIZE..(n - MIN_DOMAIN_SIZE + 1) {
        // Reset set
        for &r in &sorted {
            set[r] = 0;
        }
        // Assign group 1 (first part) and group 2 (second part)
        for &r in &sorted[..cut_pos] {
            set[r] = 1;
        }
        for &r in &sorted[cut_pos..] {
            set[r] = 2;
        }

        let tau = compute_tau(contacts, &set);
        if tau > best_tau {
            best_tau = tau;
            best_partition = set.clone();
        }
    }

    (best_tau, best_partition)
}

/// Build hierarchical domain decomposition from CA coordinates.
///
/// Returns a Vec<DomainNode> compatible with the .dat format.
pub fn build_domain_tree(
    ca: &[[f64; 3]],
    _segments: &[Segment],
) -> Vec<DomainNode> {
    let nres = ca.len();
    if nres == 0 {
        return vec![];
    }

    let contacts = build_contacts(ca);

    // Stack-based binary decomposition (like dostack in Fortran)
    let mut nodes: Vec<TreeNode> = Vec::new();

    // Root node: all residues
    nodes.push(TreeNode {
        residues: (0..nres).collect(),
        left: 0,
        right: 0,
        node_type: '*',
    });

    let mut stack: Vec<usize> = vec![0]; // indices into nodes

    while let Some(node_idx) = stack.pop() {
        let residues = nodes[node_idx].residues.clone();
        if residues.len() < 2 * MIN_DOMAIN_SIZE {
            continue;
        }

        let (tau, partition) = find_best_cut(&residues, &contacts, nres);

        // Only split if tau is significant
        if tau < 0.001 {
            continue;
        }

        // Create two child nodes
        let child1_residues: Vec<usize> = residues
            .iter()
            .filter(|&&r| partition[r] == 1)
            .copied()
            .collect();
        let child2_residues: Vec<usize> = residues
            .iter()
            .filter(|&&r| partition[r] == 2)
            .copied()
            .collect();

        if child1_residues.is_empty() || child2_residues.is_empty() {
            continue;
        }

        let left_idx = nodes.len();
        nodes.push(TreeNode {
            residues: child1_residues,
            left: 0,
            right: 0,
            node_type: '*',
        });

        let right_idx = nodes.len();
        nodes.push(TreeNode {
            residues: child2_residues,
            left: 0,
            right: 0,
            node_type: '*',
        });

        nodes[node_idx].left = left_idx + 1; // 1-based
        nodes[node_idx].right = right_idx + 1; // 1-based
        nodes[node_idx].node_type = '+';

        // Only recurse if children are large enough
        if nodes[left_idx].residues.len() >= 2 * MIN_DOMAIN_SIZE {
            stack.push(left_idx);
        }
        if nodes[right_idx].residues.len() >= 2 * MIN_DOMAIN_SIZE {
            stack.push(right_idx);
        }
    }

    // Convert to DomainNode format
    // The .dat domain tree uses the THIRD >>>> section format:
    // Each node has: index, type, left_child, right_child, nsize, nseg, (start,end) pairs
    // Where nseg and ranges describe contiguous residue segments belonging to the node.

    let mut result: Vec<DomainNode> = Vec::new();

    for (i, node) in nodes.iter().enumerate() {
        // Compute contiguous segments from residue list
        let mut sorted_res = node.residues.clone();
        sorted_res.sort();
        let segments = residues_to_segments(&sorted_res);

        let node_type = match node.node_type {
            '*' => NodeType::Root,
            '+' => NodeType::Split,
            _ => NodeType::Minus,
        };

        result.push(DomainNode {
            index: i + 1, // 1-based
            node_type,
            left_child: node.left,
            right_child: node.right,
            nseg: segments.len(),
            segments, // 1-based (start, end) pairs
        });
    }

    result
}

/// Convert sorted 0-based residue indices to 1-based (start, end) segment pairs.
fn residues_to_segments(sorted_res: &[usize]) -> Vec<(u32, u32)> {
    if sorted_res.is_empty() {
        return vec![];
    }

    let mut segments = Vec::new();
    let mut start = sorted_res[0];
    let mut end = sorted_res[0];

    for &r in &sorted_res[1..] {
        if r == end + 1 {
            end = r;
        } else {
            segments.push(((start + 1) as u32, (end + 1) as u32));
            start = r;
            end = r;
        }
    }
    segments.push(((start + 1) as u32, (end + 1) as u32));

    segments
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_residues_to_segments() {
        let residues = vec![0, 1, 2, 3, 10, 11, 12, 20];
        let segs = residues_to_segments(&residues);
        assert_eq!(segs, vec![(1, 4), (11, 13), (21, 21)]);
    }

    #[test]
    fn test_build_domain_tree_small() {
        // Small protein: just a root node
        let ca: Vec<[f64; 3]> = (0..30)
            .map(|i| [i as f64 * 3.8, 0.0, 0.0])
            .collect();
        let segments = vec![];
        let tree = build_domain_tree(&ca, &segments);
        assert!(!tree.is_empty());
        // Root should cover all residues
        assert_eq!(tree[0].index, 1);
    }
}
