use crate::types::{Segment, SseType};

/// Minimum helix length (residues) to keep without growing.
const MIN_HELIX_HARD: usize = 6;
/// Target helix length for growth (LENH parameter in Fortran).
const LENH: usize = 8;
/// Target strand length for growth (LENE parameter in Fortran).
const LENE: usize = 6;

/// Reduce 8-state DSSP to 3-state (H/E/L).
///
/// Maps: G,I,H → H; E → E; everything else → L
/// Equivalent to Fortran `simp()` with fromstring='GIHE*' tostring='HHHEL'.
fn reduce_to_3state(dssp8: &[u8]) -> Vec<u8> {
    dssp8
        .iter()
        .map(|&c| match c {
            b'G' | b'I' | b'H' => b'H',
            b'E' => b'E',
            _ => b'L',
        })
        .collect()
}

/// Grow the left border of a segment into adjacent loop residues.
///
/// Returns Ok(new_start) or Err if growth is blocked.
fn grow_left(k: usize, struc: &mut [u8]) -> Result<usize, ()> {
    if k == 0 {
        return Err(());
    }
    if struc[k - 1] != b'L' {
        return Err(());
    }
    struc[k - 1] = b' '; // mark as consumed
    Ok(k - 1)
}

/// Grow the right border of a segment into adjacent loop residues.
///
/// Returns Ok(new_end) or Err if growth is blocked.
fn grow_right(l: usize, nres: usize, struc: &mut [u8]) -> Result<usize, ()> {
    if l >= nres - 1 {
        return Err(());
    }
    if struc[l + 1] != b'L' {
        return Err(());
    }
    struc[l + 1] = b' '; // mark as consumed
    Ok(l + 1)
}

/// Assign secondary structure segments from 8-state DSSP.
///
/// Port of `puutos.f:getsecstr(nres, struc, nseg, secstr, range, lene=6, lenh=8)`.
///
/// Returns (segments, secstr_types, na, nb).
pub fn assign_segments(dssp8: &[u8]) -> (Vec<Segment>, Vec<SseType>, usize, usize) {
    let nres = dssp8.len();
    if nres == 0 {
        return (vec![], vec![], 0, 0);
    }

    // Step 1: reduce to 3-state
    let mut struc = reduce_to_3state(dssp8);

    // Step 2: find runs of same character (1-based in Fortran, 0-based here)
    let mut runs: Vec<(usize, usize, u8)> = Vec::new(); // (start, end, char) 0-based
    let mut run_start = 0;
    let mut run_char = struc[0];
    for i in 1..nres {
        if struc[i] != run_char {
            runs.push((run_start, i - 1, run_char));
            run_start = i;
            run_char = struc[i];
        }
    }
    runs.push((run_start, nres - 1, run_char));

    // Step 3: delete loop segments
    let runs: Vec<(usize, usize, u8)> = runs.into_iter().filter(|r| r.2 != b'L').collect();

    // Step 4: filter by length with growth
    let mut segments: Vec<Segment> = Vec::new();
    let mut secstr_types: Vec<SseType> = Vec::new();

    for (start, end, c) in runs {
        let mut k = start;
        let mut l = end;
        let len = l - k + 1;

        if c == b'H' && len < MIN_HELIX_HARD {
            // Too short to even try growing
            continue;
        } else if c == b'H' && len < LENH {
            // Try to grow short helix — Fortran calls growleft then growrite each iteration
            // ierr is overwritten by growrite, so only growrite failure stops the loop
            let mut ierr = false;
            while l - k + 1 < LENH && !ierr {
                // growleft — try but don't stop on failure
                if let Ok(new_k) = grow_left(k, &mut struc) {
                    k = new_k;
                }
                // growrite — failure stops the loop
                match grow_right(l, nres, &mut struc) {
                    Ok(new_l) => l = new_l,
                    Err(()) => ierr = true,
                }
            }
            if ierr && l - k + 1 < LENH {
                continue;
            }
        } else if c == b'E' && len < LENE {
            // Try to grow short strand — same pattern
            let mut ierr = false;
            while l - k + 1 < LENE && !ierr {
                if let Ok(new_k) = grow_left(k, &mut struc) {
                    k = new_k;
                }
                match grow_right(l, nres, &mut struc) {
                    Ok(new_l) => l = new_l,
                    Err(()) => ierr = true,
                }
            }
            if ierr && l - k + 1 < LENE {
                continue;
            }
        }

        // Convert to 1-based for Segment
        let sse_type = if c == b'H' {
            SseType::Helix
        } else {
            SseType::Strand
        };
        secstr_types.push(sse_type);
        segments.push(Segment {
            start: (k + 1) as u32,
            end: (l + 1) as u32,
            sse_type,
            check_start: 0,
            check_end: 0,
            checkx: 0,
        });
    }

    // Step 5: compute check ranges
    compute_check_ranges(nres, &mut segments, &secstr_types);

    let na = secstr_types.iter().filter(|&&t| t == SseType::Helix).count();
    let nb = secstr_types.iter().filter(|&&t| t == SseType::Strand).count();

    (segments, secstr_types, na, nb)
}

/// Compute check ranges for each segment.
///
/// Port of `puutos.f:getcheckrange()`.
/// Helix: x=0 if len>=9, x=1 if len 7-8, x=2 if len 6
/// Strand: x=0 if len>=5, x=1 if len 3-4, x=2 if len 2
fn compute_check_ranges(nres: usize, segments: &mut [Segment], secstr: &[SseType]) {
    for (i, seg) in segments.iter_mut().enumerate() {
        let len = (seg.end - seg.start + 1) as usize;
        let x = if secstr[i] == SseType::Helix {
            if len >= 9 {
                0
            } else if len >= 7 {
                1
            } else {
                2
            }
        } else {
            // Strand
            if len >= 5 {
                0
            } else if len >= 3 {
                1
            } else {
                2
            }
        };
        seg.checkx = x;
        seg.check_start = (seg.start as i32 - x).max(1) as u32;
        seg.check_end = (seg.end as i32 + x).min(nres as i32) as u32;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_3state_reduction() {
        // H(3)+G(3)+I(3)+E(4)+' '(1)+T(2)+S(2)+B(2)+L(2) = 22 chars
        let dssp8 = b"HHHGGGIIIEEEE TTSSBBLL";
        let result = reduce_to_3state(dssp8);
        //           HHHHHHHHHEEEE LLLLLLLL  → 3-state
        //           H(9)    E(4)L(9)        = 22
        assert_eq!(std::str::from_utf8(&result).unwrap(),
                   "HHHHHHHHHEEEELLLLLLLLL");
    }

    #[test]
    fn test_assign_segments_basic() {
        // A simple case: 10 H residues, then 5 loop, then 8 E residues
        let dssp = b"HHHHHHHHHHLLLLLEEEEEEEE";
        let (segments, secstr, na, nb) = assign_segments(dssp);
        assert_eq!(na, 1);
        assert_eq!(nb, 1);
        assert_eq!(segments.len(), 2);
        assert_eq!(secstr[0], SseType::Helix);
        assert_eq!(secstr[1], SseType::Strand);
        assert_eq!(segments[0].start, 1);
        assert_eq!(segments[0].end, 10);
        assert_eq!(segments[1].start, 16);
        assert_eq!(segments[1].end, 23);
    }

    #[test]
    fn test_short_helix_excluded() {
        // 5-residue helix should be excluded (< MIN_HELIX_HARD=6)
        let dssp = b"LLLLLHHHHHLLLLL";
        let (segments, _, na, _) = assign_segments(dssp);
        assert_eq!(na, 0);
        assert_eq!(segments.len(), 0);
    }
}
