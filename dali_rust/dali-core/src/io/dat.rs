use std::path::Path;
use ndarray::Array2;

use crate::types::{DomainNode, NodeType, Protein, Segment, SseType};

/// Parse a DaliLite .dat file into a Protein.
///
/// .dat format (Fortran fixed-width):
///   Header: `>>>> CODE  NRES NSEG  NA  NB  SECSTR`
///     - chars 5-9: code, 10-14: nres, 15-19: nseg, 20-24: na, 25-29: nb
///     - chars 32+: secondary structure characters (H/E)
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
            // First >>>>: read header, segments, CA coords
            code = line[5..10].trim().to_string();
            nres = parse_int(&line[10..15])?;
            nseg = parse_int(&line[15..20])?;
            na = parse_int(&line[20..25])?;
            nb = parse_int(&line[25..30])?;

            // Secondary structure characters starting at position 32
            if line.len() > 32 {
                let secstr_str = &line[32..line.len().min(32 + nseg)];
                secstr_chars = secstr_str.chars().collect();
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
            let ndom_header: usize = parse_int(&line[10..15])?;
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

                // Parse child pointers from the 13x area (pos 6-18)
                // These are 3 integers in i4 format: left_child, right_child, nres_or_nseg
                let child_str = if dline.len() > 6 {
                    &dline[6..dline.len().min(19)]
                } else {
                    ""
                };
                let child_ints = parse_fixed_ints_lenient(child_str, 4);
                let left_child_val = child_ints.first().copied().unwrap_or(0) as usize;
                let right_child_val = child_ints.get(1).copied().unwrap_or(0) as usize;

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
    })
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
