use std::io::Write;
use std::path::Path;
use ndarray::Array2;

use crate::types::{DomainNode, NodeType, Protein, Segment, SseType};

/// Parse a DaliLite .dat file into a Protein.
///
/// .dat format:
///   Header: `>>>> CODE NRES NSEG NA NB SECSTR`
///     - Whitespace-delimited tokens (supports codes of any length).
///     - Also handles legacy Fortran fixed-width format (5-char code at positions 5-9).
///   Segments: (6i10) — index, start, end, check_start, check_end, checkx
///   CA coords: (10f8.1) — x,y,z for each residue sequentially
///   Second `>>>>`: SSE hierarchy (skipped)
///   Third `>>>>`: domain decomposition tree
pub fn read_dat<P: AsRef<Path>>(filepath: P) -> Result<Protein, DatError> {
    let content = std::fs::read_to_string(filepath.as_ref())
        .map_err(|e| DatError::Io(e.to_string()))?;
    let lines: Vec<&str> = content.lines().collect();

    let mut idx = 0;
    let mut header_count = 0;

    let mut code = String::new();
    let mut nres = 0usize;
    let mut nseg = 0usize;
    let mut na = 0usize;
    let mut nb = 0usize;
    let mut secstr_chars: Vec<char> = Vec::new();
    let mut segments: Vec<Segment> = Vec::new();
    let mut ca = Array2::<f64>::zeros((3, 0));
    let mut domain_tree: Vec<DomainNode> = Vec::new();
    let mut sequence = String::new();

    while idx < lines.len() {
        let line = lines[idx];

        if !line.starts_with(">>>>") {
            // Check for -sequence line
            if line.starts_with("-sequence") {
                // Format: -sequence "SEQUENCE..."
                if let Some(start) = line.find('"') {
                    if let Some(end) = line.rfind('"') {
                        if end > start {
                            sequence = line[start + 1..end].to_string();
                        }
                    }
                }
            }
            idx += 1;
            continue;
        }

        header_count += 1;

        if header_count == 1 {
            // First >>>>: read header, segments, CA coords.
            // Whitespace-delimited: >>>> CODE NRES NSEG NA NB SECSTR
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 6 {
                code = parts[1].to_string();
                nres = parts[2].parse().map_err(|_| DatError::Format(format!("Bad nres '{}'", parts[2])))?;
                nseg = parts[3].parse().map_err(|_| DatError::Format(format!("Bad nseg '{}'", parts[3])))?;
                na = parts[4].parse().map_err(|_| DatError::Format(format!("Bad na '{}'", parts[4])))?;
                nb = parts[5].parse().map_err(|_| DatError::Format(format!("Bad nb '{}'", parts[5])))?;
                if parts.len() > 6 {
                    secstr_chars = parts[6].chars().take(nseg).collect();
                }
            } else {
                return Err(DatError::Format(format!("Bad >>>> header: '{}'", line)));
            }

            idx += 1;

            // Read nseg segment lines: format (6i10)
            segments = Vec::with_capacity(nseg);
            for _ in 0..nseg {
                if idx >= lines.len() {
                    return Err(DatError::Format("Unexpected end of file in segments".into()));
                }
                let seg_line = lines[idx];
                idx += 1;

                let vals = parse_fixed_ints(seg_line, 10, 6)?;
                let sse_idx = vals[0] as usize;
                let sse_type = if sse_idx <= secstr_chars.len() {
                    SseType::from_char(secstr_chars[sse_idx - 1]).unwrap_or(SseType::Helix)
                } else {
                    SseType::Helix
                };

                segments.push(Segment {
                    start: vals[1] as u32,
                    end: vals[2] as u32,
                    sse_type,
                    check_start: vals[3] as u32,
                    check_end: vals[4] as u32,
                    checkx: vals[5] as i32,
                });
            }

            // Read CA coordinates: format (10f8.1)
            let total_vals = 3 * nres;
            let mut ca_vals: Vec<f64> = Vec::with_capacity(total_vals);
            while ca_vals.len() < total_vals && idx < lines.len() {
                let coord_line = lines[idx];
                idx += 1;
                let remaining = total_vals - ca_vals.len();
                let nfloats = remaining.min(10);
                for k in 0..nfloats {
                    let start = k * 8;
                    let end = (start + 8).min(coord_line.len());
                    if start >= coord_line.len() {
                        break;
                    }
                    let s = coord_line[start..end].trim();
                    let val: f64 = s.parse().map_err(|_| {
                        DatError::Format(format!("Bad float '{}' in CA coords", s))
                    })?;
                    ca_vals.push(val);
                }
            }

            if ca_vals.len() != total_vals {
                return Err(DatError::Format(format!(
                    "Expected {} CA values, got {}",
                    total_vals,
                    ca_vals.len()
                )));
            }

            // Reshape: file order is x1,y1,z1,x2,y2,z2,... -> (3, nres)
            ca = Array2::zeros((3, nres));
            for j in 0..nres {
                ca[[0, j]] = ca_vals[j * 3];
                ca[[1, j]] = ca_vals[j * 3 + 1];
                ca[[2, j]] = ca_vals[j * 3 + 2];
            }
            continue;
        } else if header_count == 2 {
            // Second >>>>: SSE hierarchy (skip)
            idx += 1;
            continue;
        } else if header_count == 3 {
            // Third >>>>: domain decomposition tree
            // Whitespace: >>>> CODE NDOM  (or legacy fixed: >>>> at 0-4, code at 5-9, ndom at 10-14)
            let parts: Vec<&str> = line.split_whitespace().collect();
            let ndom_header: usize = if parts.len() >= 3 {
                parts[2].parse().unwrap_or(0)
            } else if line.len() >= 15 {
                parse_int(&line[10..15]).unwrap_or(0)
            } else {
                0
            };
            idx += 1;

            domain_tree = Vec::with_capacity(ndom_header);
            for _ in 0..ndom_header {
                if idx >= lines.len() {
                    break;
                }
                let dline = lines[idx];
                idx += 1;

                if dline.len() < 6 {
                    continue;
                }

                // Parse: i4 (pos 0-3), 1x (pos 4), a1 (pos 5), 13x (pos 6-18)
                let j: usize = dline[0..4].trim().parse().map_err(|_| {
                    DatError::Format(format!("Bad domain index in '{}'", &dline[0..4]))
                })?;
                let node_char = dline.chars().nth(5).unwrap_or('-');
                let node_type = NodeType::from_char(node_char);

                // Child pointers: Fortran (i4,1x,a1,1x,3i4) puts them at pos 7-10 and
                // 11-14 (pos 15-18 is the node size). Files from the pre-2026-10 Rust
                // writer have them one column left (6-9, 10-13); trimming the 7-10 /
                // 11-14 windows reads those correctly too for children < 1000.
                let field = |a: usize, b: usize| -> usize {
                    dline.get(a..b.min(dline.len()))
                        .and_then(|s| s.trim().parse().ok())
                        .unwrap_or(0)
                };
                let left_child_val = field(7, 11);
                let right_child_val = field(11, 15);

                // 400i4 starting at position 19: domns + segment ranges
                // ints[0] = nseg (number of segments in this domain)
                // ints[1:] = (start, end) pairs for each segment
                let int_str = if dline.len() > 19 { &dline[19..] } else { "" };
                let ints = parse_fixed_ints_lenient(int_str, 4);

                let (left_child, right_child, dom_nseg, dom_segs) = if !ints.is_empty() {
                    let ns = ints[0] as usize;
                    let mut segs = Vec::new();
                    for s in 0..ns {
                        if 1 + s * 2 + 1 < ints.len() {
                            segs.push((ints[1 + s * 2] as u32, ints[1 + s * 2 + 1] as u32));
                        }
                    }
                    (left_child_val, right_child_val, ns, segs)
                } else {
                    (0, 0, 0, Vec::new())
                };

                domain_tree.push(DomainNode {
                    index: j,
                    node_type,
                    left_child,
                    right_child,
                    nseg: dom_nseg,
                    segments: dom_segs,
                });
            }
            break;
        }

        idx += 1;
    }

    let secstr: Vec<SseType> = secstr_chars
        .iter()
        .map(|&c| SseType::from_char(c).unwrap_or(SseType::Helix))
        .collect();

    // .dat files don't store PDB residue IDs; use sequential 1..=nres
    let resid_map: Vec<i32> = (1..=nres as i32).collect();

    Ok(Protein {
        code,
        nres,
        nseg,
        na,
        nb,
        segments,
        secstr,
        ca,
        sequence,
        domain_tree,
        resid_map,
        numbering: crate::ResidNumbering::Sequential,
    })
}

/// One node of the SSE-level domain tree (the second `>>>>` section of a .dat file).
///
/// This is the tree PARSI reads: each node lists the SSEs (1-based segment indices,
/// ascending) that belong to it.
#[derive(Debug, Clone, PartialEq)]
pub struct SseTreeNode {
    pub node_type: char,
    pub left_child: usize,  // 1-based, 0 = none
    pub right_child: usize, // 1-based, 0 = none
    pub sses: Vec<usize>,
}

/// Derive the SSE-level domain tree from the residue-level tree.
///
/// Port of puutos.f `assignsegments` + `prunetree` + `getdomseglist`:
///   1. each residue belongs to the deepest (highest-index) node covering it
///   2. each SSE goes to the node owning its midpoint residue
///   3. parents collect their children's SSEs
///   4. empty nodes are deleted (the parent absorbs the sibling's children)
///   5. surviving nodes are renumbered top-down, preserving order
pub fn derive_sse_tree(
    nres: usize,
    segments: &[Segment],
    domain_tree: &[DomainNode],
) -> Vec<SseTreeNode> {
    let ndom = domain_tree.len();
    if ndom == 0 {
        return Vec::new();
    }

    // 1-based node arrays, indexed by position in domain_tree (+1)
    let mut child = vec![[0usize; 2]; ndom + 1];
    let mut node_type = vec![' '; ndom + 1];
    for (i, node) in domain_tree.iter().enumerate() {
        let idom = i + 1;
        let valid = |c: usize| if c <= ndom { c } else { 0 };
        child[idom] = [valid(node.left_child), valid(node.right_child)];
        node_type[idom] = node.node_type.to_char();
    }
    let mut parent = vec![0usize; ndom + 1];
    for idom in 1..=ndom {
        for &c in &child[idom] {
            if c > 0 {
                parent[c] = idom;
            }
        }
    }

    // assignsegments: bottom node per residue, then segment midpoints
    let mut bottom = vec![0usize; nres + 1];
    for (i, node) in domain_tree.iter().enumerate() {
        for &(s, e) in &node.segments {
            for ires in (s as usize).max(1)..=(e as usize).min(nres) {
                bottom[ires] = i + 1;
            }
        }
    }
    let mut node_sses: Vec<Vec<usize>> = vec![Vec::new(); ndom + 1];
    for (iseg, seg) in segments.iter().enumerate() {
        let mid = ((seg.start + seg.end) / 2) as usize;
        let idom = if mid <= nres { bottom[mid] } else { 0 };
        if idom > 0 {
            node_sses[idom].push(iseg + 1);
        }
    }

    // prunetree: build parent ranges bottom-up
    for idom in (1..=ndom).rev() {
        let pnode = parent[idom];
        if pnode > 0 && child[pnode][0] == idom {
            let [inode, jnode] = child[pnode];
            let mut add = node_sses[inode].clone();
            if jnode > 0 {
                add.extend_from_slice(&node_sses[jnode]);
            }
            node_sses[pnode].extend(add);
        }
    }

    // delete empty nodes
    let mut active = vec![true; ndom + 1];
    for idom in (1..=ndom).rev() {
        let pnode = parent[idom];
        if !(active[idom] && pnode > 0) {
            continue;
        }
        let [c1, c2] = child[pnode];
        let jnode = if c2 == idom { c1 } else { c2 }; // sibling
        if node_sses[idom].is_empty() && jnode > 0 {
            child[pnode] = child[jnode];
            active[jnode] = false;
        }
    }
    for idom in 1..=ndom {
        if !active[idom] || node_sses[idom].is_empty() {
            parent[idom] = 0;
            child[idom] = [0, 0];
        }
    }
    for idom in 1..=ndom {
        for &c in &child[idom].clone() {
            if c > 0 {
                parent[c] = idom;
            }
        }
    }

    // minimal tree top-down
    let mut list = vec![1usize];
    let mut tsil = vec![0usize; ndom + 1];
    tsil[1] = 1;
    for idom in 2..=ndom {
        if parent[idom] > 0 {
            list.push(idom);
            tsil[idom] = list.len();
        }
    }

    list.iter()
        .map(|&idom| {
            let mut sses = node_sses[idom].clone();
            sses.sort_unstable();
            sses.dedup();
            SseTreeNode {
                node_type: node_type[idom],
                left_child: tsil[child[idom][0]],
                right_child: tsil[child[idom][1]],
                sses,
            }
        })
        .collect()
}

/// Section header: Fortran `('>>>> ',a4,a1,...)` layout when the code fits in 5
/// characters, whitespace-delimited otherwise (both readers tokenize on whitespace).
fn section_header(code: &str, fields: &[usize], tail: &str) -> String {
    let fixed = code.len() <= 5
        && !code.contains(char::is_whitespace)
        && fields.iter().all(|&v| v < 10_000);
    let mut s = String::from(">>>> ");
    if fixed {
        s.push_str(&format!("{:<5}", code));
        for v in fields {
            s.push_str(&format!("{:>5}", v));
        }
        if !tail.is_empty() {
            s.push_str("  ");
            s.push_str(tail);
        }
    } else {
        s.push_str(code);
        for v in fields {
            s.push_str(&format!(" {}", v));
        }
        if !tail.is_empty() {
            s.push(' ');
            s.push_str(tail);
        }
    }
    s
}

/// Write a Protein to a DaliLite .dat file.
///
/// Produces the same 3 `>>>>` sections as Fortran puutos.f:
///   Section 1: header + segments (6i10) + CA coords (10f8.1)
///   Section 2: SSE-level domain tree, format (i4,1x,a1,1x,3i4,200i4) — read by PARSI
///   Section 3: residue-level domain tree, same format with node size + ranges
pub fn write_dat<P: AsRef<Path>>(protein: &Protein, filepath: P) -> Result<(), DatError> {
    let mut f = std::fs::File::create(filepath.as_ref())
        .map_err(|e| DatError::Io(e.to_string()))?;
    let io = |e: std::io::Error| DatError::Io(e.to_string());

    let code = &protein.code;
    let secstr_str: String = protein.secstr.iter().map(|s| s.to_char()).collect();

    // Section 1: header, segments, CA
    writeln!(f, "{}", section_header(
        code, &[protein.nres, protein.nseg, protein.na, protein.nb], &secstr_str)).map_err(io)?;

    for (i, seg) in protein.segments.iter().enumerate() {
        writeln!(f, "{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}",
                 i + 1, seg.start, seg.end, seg.check_start, seg.check_end, seg.checkx)
            .map_err(io)?;
    }

    // CA coordinates: 10f8.1 format (x1,y1,z1,x2,y2,z2,...)
    let total_vals = 3 * protein.nres;
    let mut val_idx = 0;
    while val_idx < total_vals {
        let nfloats = (total_vals - val_idx).min(10);
        let mut line = String::with_capacity(80);
        for k in 0..nfloats {
            let res = (val_idx + k) / 3;
            let coord = (val_idx + k) % 3;
            line.push_str(&format!("{:>8.1}", protein.ca[[coord, res]]));
        }
        writeln!(f, "{}", line).map_err(io)?;
        val_idx += nfloats;
    }

    // Section 2: SSE-level domain tree
    let sse_tree = derive_sse_tree(protein.nres, &protein.segments, &protein.domain_tree);
    writeln!(f, "{}", section_header(code, &[sse_tree.len()], "")).map_err(io)?;
    for (i, node) in sse_tree.iter().enumerate() {
        let mut line = format!("{:>4} {} {:>4}{:>4}{:>4}",
            i + 1, node.node_type, node.left_child, node.right_child, node.sses.len());
        for &s in &node.sses {
            line.push_str(&format!("{:>4}", s));
        }
        writeln!(f, "{}", line).map_err(io)?;
    }

    // Section 3: residue-level domain tree
    writeln!(f, "{}", section_header(code, &[protein.domain_tree.len()], "")).map_err(io)?;
    for node in &protein.domain_tree {
        let size: u32 = node.segments.iter().map(|&(s, e)| e + 1 - s).sum();
        let mut line = format!("{:>4} {} {:>4}{:>4}{:>4}{:>4}",
            node.index, node.node_type.to_char(), node.left_child, node.right_child,
            size, node.segments.len());
        for &(start, end) in &node.segments {
            line.push_str(&format!("{:>4}{:>4}", start, end));
        }
        writeln!(f, "{}", line).map_err(io)?;
    }

    Ok(())
}

/// Parse a fixed-width integer field.
fn parse_int(s: &str) -> Result<usize, DatError> {
    s.trim()
        .parse()
        .map_err(|_| DatError::Format(format!("Bad integer '{}'", s.trim())))
}

/// Parse fixed-width integers from a line.
fn parse_fixed_ints(line: &str, width: usize, count: usize) -> Result<Vec<i64>, DatError> {
    let mut vals = Vec::with_capacity(count);
    for k in 0..count {
        let start = k * width;
        let end = (start + width).min(line.len());
        if start >= line.len() {
            return Err(DatError::Format(format!(
                "Line too short for {} ints of width {}: '{}'",
                count, width, line
            )));
        }
        let s = line[start..end].trim();
        let val: i64 = s.parse().map_err(|_| {
            DatError::Format(format!("Bad integer '{}' in '{}'", s, line))
        })?;
        vals.push(val);
    }
    Ok(vals)
}

/// Parse fixed-width integers, stopping at parse failure (lenient).
fn parse_fixed_ints_lenient(line: &str, width: usize) -> Vec<i64> {
    let mut vals = Vec::new();
    let mut pos = 0;
    let trimmed = line.trim_end_matches('\n');
    while pos + width <= trimmed.len() {
        let s = trimmed[pos..pos + width].trim();
        if s.is_empty() {
            pos += width;
            continue;
        }
        match s.parse::<i64>() {
            Ok(v) => vals.push(v),
            Err(_) => break,
        }
        pos += width;
    }
    vals
}

#[derive(Debug)]
pub enum DatError {
    Io(String),
    Format(String),
}

impl std::fmt::Display for DatError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DatError::Io(msg) => write!(f, "I/O error: {}", msg),
            DatError::Format(msg) => write!(f, "Format error: {}", msg),
        }
    }
}

impl std::error::Error for DatError {}
