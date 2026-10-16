//! PUU domain decomposition — port of DaliLite puu.f (L. Holm, C. Sander, Proteins
//! 19:256-268, 1994), single-chain path.
//!
//! Produces the unfolding tree that puu writes to `units.puu` and puutos reads.
//! Faithful to the Fortran, including its quirks, because DALI's PARSI walks this
//! tree and DP/FILTER95 read its node types:
//!
//! - all arithmetic in f32 (gfortran -O3, no FMA, no reassociation), sums in the
//!   Fortran loop order
//! - the spectral cut's "random" start is gfortran's RAN intrinsic reseeded on every
//!   call (ran.f is not linked), i.e. a constant vector; see `gfortran_ran`
//! - getcoor reads coordinates list-directed: a field of -100 or below runs into
//!   its neighbour and the unread coordinates keep stale values (see read_pdb_atoms)
//! - getcoor's backbone-incomplete rollback leaves the dropped residue's first atom
//!   behind, attributed to the next residue (or, at end of file, to residue nres+1,
//!   where it only affects centering)
//! - in the single-chain path contstren is passed the per-chain atom-count array,
//!   which is never filled (static, zero); every domain has natom=1 and
//!   "compactness" is 100 x the intra-unit contact count
//!
//! The multi-chain branch of doprotein (chain wildcard '*') is not ported.

use std::io::Read;
use std::path::Path;

/// One atom as getcoor keeps it.
#[derive(Debug, Clone)]
pub struct PuuAtom {
    pub name: [u8; 3], // columns 14-16
    pub res: usize,    // 1-based residue index; nres+1 for an end-of-file rollback stray
    pub xyz: [f32; 3],
}

/// getcoor output.
#[derive(Debug, Clone)]
pub struct PuuStructure {
    pub nres: usize,
    pub atoms: Vec<PuuAtom>,
    pub resno: Vec<String>, // per residue: chain + right-adjusted resseq/icode (5 chars)
}

/// The subset of a DSSP file puu reads.
#[derive(Debug, Clone)]
pub struct DsspInput {
    pub struc: Vec<u8>,       // column 17
    pub bp: Vec<[i32; 2]>,    // bridge partners, sequential 1-based, 0 = none
}

/// puu parameters (puu.f `defaults` + the puu.default mpidali.pm writes).
#[derive(Debug, Clone)]
pub struct PuuParams {
    pub breakdist: f32,
    pub minseglen: usize,
    pub hbond: f32,
    pub compact: f32,
    pub taucutoff: f32,
    pub seed: i32,
    pub lsplit: bool,
    pub printdepth: i32,
}

impl Default for PuuParams {
    fn default() -> Self {
        PuuParams {
            breakdist: 5.0,
            minseglen: 10,
            hbond: 15.0,
            compact: 80.0,
            taucutoff: 2.6,
            seed: 1_234_567_811,
            lsplit: true,     // mpidali.pm puu.default
            printdepth: -99,  // mpidali.pm puu.default
        }
    }
}

/// One node of the unfolding tree, as written to units.puu.
#[derive(Debug, Clone, PartialEq)]
pub struct PuuNode {
    pub child: [usize; 2], // 1-based, 0 = none
    pub node_type: char,   // '+', '*', '=', '-'
    pub tau: f32,
    pub tert: f32,
    pub nbp: i32,
    pub size: usize,
    pub ranges: Vec<(usize, usize)>, // 1-based residue ranges
}

#[derive(Debug)]
pub enum PuuError {
    Io(String),
    Empty,
}

impl std::fmt::Display for PuuError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PuuError::Io(s) => write!(f, "{}", s),
            PuuError::Empty => write!(f, "no residues"),
        }
    }
}

impl std::error::Error for PuuError {}

const MAXGRID: i32 = 20;
const MAXNSEG: usize = 25;

// ---------------------------------------------------------------------------
// Input: getcoor (util.f) and getdssp (util.f)
// ---------------------------------------------------------------------------

fn read_text<P: AsRef<Path>>(path: P) -> Result<String, PuuError> {
    let p = path.as_ref();
    let file = std::fs::File::open(p).map_err(|e| PuuError::Io(format!("{}: {}", p.display(), e)))?;
    let mut s = String::new();
    let r: Result<usize, std::io::Error> = if p.extension().map_or(false, |x| x == "gz") {
        flate2::read::GzDecoder::new(file).read_to_string(&mut s)
    } else {
        let mut f = file;
        f.read_to_string(&mut s)
    };
    r.map_err(|e| PuuError::Io(format!("{}: {}", p.display(), e)))?;
    Ok(s)
}

/// Fortran `a80`/`a136` read: truncate or pad with blanks, byte-indexed.
fn fixed(line: &str, width: usize) -> Vec<u8> {
    let mut b: Vec<u8> = line.bytes().take(width).collect();
    b.resize(width, b' ');
    b
}

/// util.f onelettercode: '-' marks an unacceptable residue.
fn onelettercode(aaa: &[u8]) -> u8 {
    const AMINOACID: &[u8] = b"ALAARGASNASPCYSGLUGLNGLYHISILELEULYSMETPHEPROSERTHRTRPTYRVAL\
ASXGLXACDALBALIABUAROBASBETHSEHYPHYLORNPCASARTAUTHYUNKACEFOR\
CYHCSHCSSCYXILUPRZPR0CPRTRYHOHSEMMSE";
    const AASYMBOL: &[u8] = b"ARNDCEQGHILKMFPSTWYVBZXXXXXXXXXXXXXXXX--CCCCIPPPW-XM";
    let mut i = 0;
    while i + 2 < 154 {
        if &AMINOACID[i..i + 3] == aaa {
            return AASYMBOL[i / 3];
        }
        i += 3;
    }
    b'-'
}

/// Port of util.f getcoor for one chain (`chainid`), lhetatm as given.
pub fn read_pdb_atoms<P: AsRef<Path>>(
    path: P,
    chainid: u8,
    lhetatm: bool,
) -> Result<PuuStructure, PuuError> {
    let text = read_text(path)?;
    let mut atoms: Vec<PuuAtom> = Vec::new();
    let mut resno: Vec<String> = Vec::new();
    let mut nres = 0usize;
    let mut oldresname: Vec<u8> = b"!@#$  ".to_vec();
    let mut slots: Vec<[f32; 3]> = Vec::new(); // setup_protein's xyz array, by atom slot
    let (mut ln, mut lc, mut lca, mut lo) = (true, true, true, true);

    // Drop residue nres (incomplete backbone). Fortran sets natom to the index of the
    // residue's FIRST atom, so that atom survives and belongs to the next residue.
    let rollback = |atoms: &mut Vec<PuuAtom>, nres: usize| {
        let mut i = atoms.len() + 1; // Fortran 1-based i = natom+1
        while i > 1 && atoms[i - 2].res == nres {
            i -= 1;
        }
        atoms.truncate(i); // natom = i
    };

    for raw in text.lines() {
        let line = fixed(raw, 80);
        if &line[0..6] == b"ENDMDL" {
            break;
        }
        let is_atom = &line[0..6] == b"ATOM  ";
        let is_het = &line[0..6] == b"HETATM";
        if !(is_atom || (is_het && lhetatm)) {
            continue;
        }
        if chainid != b'*' && line[21] != chainid {
            continue;
        }
        let resname = line[21..27].to_vec();
        let altloc = line[16];
        if altloc != b' ' && altloc != b'A' && altloc != b'1' {
            continue;
        }
        let aaa = &line[17..20];
        if onelettercode(aaa) == b'-' && !is_het {
            continue;
        }
        let d = line[13];
        if d == b'H' || d == b'D' {
            continue;
        }
        if resname != oldresname {
            if nres > 0 && (ln || lo || lc || lca) {
                rollback(&mut atoms, nres);
                resno.pop();
                nres -= 1;
            }
            ln = true;
            lo = true;
            lc = true;
            lca = true;
            nres += 1;
            let mut rno = line[22..26].to_vec();
            let ic = line[26];
            if ic != b' ' {
                rno.remove(0);
                rno.push(ic);
            }
            let mut r = String::with_capacity(5);
            r.push(line[21] as char);
            r.push_str(&String::from_utf8_lossy(&rno));
            resno.push(r);
            oldresname = resname;
        }
        let name = [line[13], line[14], line[15]];
        match &name {
            b"N  " => ln = false,
            b"CA " => lca = false,
            b"C  " => lc = false,
            b"O  " => lo = false,
            _ => {}
        }
        // List-directed read of columns 31-54, as getcoor does. A coordinate of
        // -100 or below fills its 8-column field and runs into the previous one
        // ("48.118-100.952"); the read then fails at that token: values before it
        // are assigned, the rest keep whatever this atom slot of the (static,
        // zero-initialised) coordinate array last held. dsspcmbi's fscanf reads
        // such lines correctly, so only puu sees the stale coordinates.
        let slot = atoms.len();
        if slots.len() <= slot {
            slots.resize(slot + 1, [0.0; 3]);
        }
        let field = String::from_utf8_lossy(&line[30..54]).to_string();
        let mut toks = field.split(|c: char| c.is_whitespace() || c == ',').filter(|t| !t.is_empty());
        for k in 0..3 {
            match toks.next().and_then(|t| t.parse::<f32>().ok()) {
                Some(v) => slots[slot][k] = v,
                None => break,
            }
        }
        atoms.push(PuuAtom { name, res: nres, xyz: slots[slot] });
    }
    if nres > 0 && (ln || lo || lc || lca) {
        rollback(&mut atoms, nres);
        resno.pop();
        nres -= 1;
    }
    if nres == 0 {
        return Err(PuuError::Empty);
    }
    Ok(PuuStructure { nres, atoms, resno })
}

/// Port of util.f getdssp for what puu uses: secondary structure and bridge partners
/// renumbered to sequential residue indices (-99 = partner in another chain).
pub fn read_dssp<P: AsRef<Path>>(path: P, chainid: u8) -> Result<DsspInput, PuuError> {
    let text = read_text(path)?;
    let mut lines = text.lines();
    let mut ntot = 99_999usize;
    for raw in lines.by_ref() {
        let line = fixed(raw, 136);
        if &line[18..42] == b"TOTAL NUMBER OF RESIDUES" {
            if let Some(n) = String::from_utf8_lossy(&line).split_whitespace().next()
                .and_then(|s| s.parse().ok()) {
                ntot = n;
            }
        }
        if &line[0..4] == b"  # " {
            break;
        }
    }
    let mut struc = Vec::new();
    let mut raw_bp: Vec<[i32; 2]> = Vec::new();
    let mut fsspresno: std::collections::HashMap<i32, i32> = std::collections::HashMap::new();
    let mut linecount = 0;
    for raw in lines {
        let line = fixed(raw, 136);
        if line[13] != b'!' {
            linecount += 1;
        }
        if linecount > ntot {
            break;
        }
        if line[13] == b'!' {
            continue;
        }
        let i: i32 = String::from_utf8_lossy(&line[0..5]).trim().parse().unwrap_or(0);
        if chainid != b'*' && line[11] != chainid {
            fsspresno.insert(i, 0);
            continue;
        }
        struc.push(line[16]);
        let p = |a: usize, b: usize| String::from_utf8_lossy(&line[a..b]).trim().parse::<i32>().unwrap_or(0);
        raw_bp.push([p(25, 29), p(29, 33)]);
        fsspresno.insert(i, struc.len() as i32);
    }
    let bp = raw_bp
        .iter()
        .map(|b| {
            let mut out = [0i32; 2];
            for k in 0..2 {
                out[k] = if b[k] > 0 {
                    match fsspresno.get(&b[k]).copied().unwrap_or(0) {
                        x if x > 0 => x,
                        _ => -99,
                    }
                } else {
                    b[k]
                };
            }
            out
        })
        .collect();
    Ok(DsspInput { struc, bp })
}

// ---------------------------------------------------------------------------
// Contacts: getcontact_shortlist, center, fillgrid, init_search, addhbond
// ---------------------------------------------------------------------------

#[inline]
fn nint(x: f32) -> i32 {
    x.round() as i32
}

#[inline]
fn fung(x: f32) -> i32 {
    nint(x / 2.0)
}

#[inline]
fn distance(a: &[f32; 3], b: &[f32; 3]) -> f32 {
    ((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2])).sqrt()
}

#[inline]
fn lgrid(gx: i32, gy: i32, gz: i32) -> bool {
    gx.abs() <= MAXGRID && gy.abs() <= MAXGRID && gz.abs() <= MAXGRID
}

#[inline]
fn contact_strength(r: f32) -> f32 {
    if r <= 4.0 { 1.0 } else { 0.0 }
}

fn center(xyz: &mut [[f32; 3]]) {
    if xyz.is_empty() {
        return;
    }
    let mut c = [0.0f32; 3];
    for a in xyz.iter() {
        for j in 0..3 {
            c[j] += a[j];
        }
    }
    let n = xyz.len() as f32;
    for j in 0..3 {
        c[j] /= n;
    }
    for a in xyz.iter_mut() {
        for j in 0..3 {
            a[j] -= c[j];
        }
    }
}

/// init_search(4.0): grid-cell offsets to scan per atom, in Fortran order.
fn init_search(realcutoff: f32) -> Vec<[i32; 3]> {
    let cutoff = nint(realcutoff + 0.5);
    let cutoff2 = nint(realcutoff * realcutoff);
    let r = cutoff + 2;
    let mut test = std::collections::HashSet::new();
    let (mut mn, mut mx) = ([0i32; 3], [0i32; 3]);
    for dx in -r..=r {
        for dy in -r..=r {
            for dz in -r..=r {
                let g = [fung(dx as f32), fung(dy as f32), fung(dz as f32)];
                for k in 0..3 {
                    mn[k] = mn[k].min(g[k]);
                    mx[k] = mx[k].max(g[k]);
                }
                let m2 = |g: i32| {
                    let g = g * 2;
                    ((g - 2) * (g - 2)).min((g + 2) * (g + 2)).min(g * g)
                };
                if m2(g[0]) + m2(g[1]) + m2(g[2]) < cutoff2 {
                    test.insert(g);
                }
            }
        }
    }
    let mut out = Vec::new();
    for dx in mn[0]..=mx[0] {
        for dy in mn[1]..=mx[1] {
            for dz in mn[2]..=mx[2] {
                if test.contains(&[dx, dy, dz]) {
                    out.push([dx, dy, dz]);
                }
            }
        }
    }
    out
}

/// Shortlisted contacts (point list order matters: it drives the f32 power iteration).
struct Contacts {
    point: Vec<[usize; 2]>, // 1-based residues
    occx: Vec<f32>,
}

fn get_contacts(s: &PuuStructure) -> Contacts {
    let nres = s.nres;
    let natom = s.atoms.len();
    let mut xyz: Vec<[f32; 3]> = s.atoms.iter().map(|a| a.xyz).collect();
    let res: Vec<usize> = s.atoms.iter().map(|a| a.res).collect();

    center(&mut xyz);

    // residues with any atom outside the grid go to "set 2"
    let mut lskip = vec![false; nres + 2];
    for i in 0..natom {
        let r = res[i];
        if !lskip[r] && !lgrid(fung(xyz[i][0]), fung(xyz[i][1]), fung(xyz[i][2])) {
            lskip[r] = true;
        }
    }
    let lset2 = (1..=nres).any(|r| lskip[r]);
    let active: Vec<bool> = (0..natom).map(|i| !lskip[res[i]]).collect();

    // fillgrid centers again before binning; distances use these coordinates
    center(&mut xyz);
    let side = (2 * MAXGRID + 1) as usize;
    let cell = |g: [i32; 3]| -> usize {
        (((g[0] + MAXGRID) as usize * side) + (g[1] + MAXGRID) as usize) * side + (g[2] + MAXGRID) as usize
    };
    let mut grid: Vec<Vec<usize>> = vec![Vec::new(); side * side * side];
    for i in 0..natom {
        if active[i] {
            let g = [fung(xyz[i][0]), fung(xyz[i][1]), fung(xyz[i][2])];
            if lgrid(g[0], g[1], g[2]) {
                grid[cell(g)].push(i);
            }
        }
    }
    let search = init_search(4.0);

    let mut point = Vec::new();
    let mut occx = Vec::new();
    let save = |r: usize, cont: &[f32], point: &mut Vec<[usize; 2]>, occx: &mut Vec<f32>| {
        if r == 0 || r > nres {
            return; // end-of-file stray residue: its row is never used downstream
        }
        for j in 1..=nres {
            if cont[j] > 0.0 {
                point.push([r, j]);
                occx.push(cont[j]);
            }
        }
    };

    let mut cont = vec![0.0f32; nres + 2];
    let mut oldres = 0usize;
    for i in 0..natom {
        if !active[i] {
            continue;
        }
        let ires = res[i];
        if ires != oldres {
            if oldres > 0 && !lskip[oldres] {
                save(oldres, &cont, &mut point, &mut occx);
            }
            cont.iter_mut().for_each(|c| *c = 0.0);
            oldres = ires;
        }
        let g = [fung(xyz[i][0]), fung(xyz[i][1]), fung(xyz[i][2])];
        for off in &search {
            let g1 = [g[0] + off[0], g[1] + off[1], g[2] + off[2]];
            if !lgrid(g1[0], g1[1], g1[2]) {
                continue;
            }
            for &j in &grid[cell(g1)] {
                cont[res[j]] += contact_strength(distance(&xyz[i], &xyz[j]));
            }
        }
        if lset2 {
            for j in 0..natom {
                if lskip[res[j]] {
                    cont[res[j]] += contact_strength(distance(&xyz[i], &xyz[j]));
                }
            }
        }
    }
    if oldres > 0 && !lskip[oldres] {
        save(oldres, &cont, &mut point, &mut occx);
    }

    if lset2 {
        let mut first = vec![usize::MAX; nres + 2];
        let mut last = vec![0usize; nres + 2];
        for i in 0..natom {
            let r = res[i];
            if first[r] == usize::MAX {
                first[r] = i;
            }
            last[r] = i;
        }
        for ires in 1..=nres {
            if !lskip[ires] {
                continue;
            }
            cont.iter_mut().for_each(|c| *c = 0.0);
            for jres in 1..=nres {
                if first[jres] == usize::MAX {
                    continue;
                }
                if distance(&xyz[first[ires]], &xyz[first[jres]]) < 25.0 {
                    for i in first[ires]..=last[ires] {
                        for j in first[jres]..=last[jres] {
                            cont[jres] += contact_strength(distance(&xyz[i], &xyz[j]));
                        }
                    }
                }
            }
            save(ires, &cont, &mut point, &mut occx);
        }
    }

    Contacts { point, occx }
}

/// setup_protein's DSSP part: bridge partners made symmetric, then H-bond weights
/// added to the contact list.
fn add_hbonds(c: &mut Contacts, dssp: &DsspInput, nres: usize, hbond: f32) -> Vec<[i32; 2]> {
    let nres1 = dssp.bp.len();
    let mut bp: Vec<[i32; 2]> = vec![[0, 0]; nres1.max(nres) + 1]; // 1-based
    for i in 1..=nres1 {
        let mut b = dssp.bp[i - 1];
        for k in 0..2 {
            if b[k] == -99 {
                b[k] = 0;
            }
        }
        bp[i] = b;
    }
    // checkdssp
    for i in 1..=nres1 {
        for k in 0..2 {
            let j = bp[i][k];
            if j > 0 && (j as usize) < bp.len() {
                let j = j as usize;
                if bp[j][0] != i as i32 && bp[j][1] != i as i32 {
                    if bp[j][0] == 0 {
                        bp[j][0] = i as i32;
                    } else if bp[j][1] == 0 {
                        bp[j][1] = i as i32;
                    }
                }
            }
        }
    }
    // addhbond
    for i in 1..=nres1 {
        for k in 0..2 {
            let j = bp[i][k];
            if j <= 0 || i > nres || j as usize > nres {
                continue;
            }
            let pj = [i, j as usize];
            match c.point.iter().position(|p| *p == pj) {
                Some(ix) => c.occx[ix] += hbond,
                None => {
                    c.point.push(pj);
                    c.occx.push(hbond);
                }
            }
        }
    }
    bp
}

// ---------------------------------------------------------------------------
// Spectral bisection: correspondence_shortlist, schmidt, definesets, tau, savecut
// ---------------------------------------------------------------------------

/// What `ran(seed)` actually is in the DaliLite build. The Makefile links puu from
/// main.f puu.f util.f -- not ran.f -- so `ran` resolves to the gfortran intrinsic
/// RAN(flag) = RAND(flag), which treats any flag other than 0/1 as a NEW SEED: every
/// call reseeds with the same value and returns the same number, and the seed variable
/// is not updated. libgfortran rand.c: Park-Miller step, then (s-1) with the low 9 bits
/// cleared, scaled by 2^-31. For the default seed this is 0.18356466.
///
/// Consequence for puu: the "random" start of the second eigenvector is constant, so
/// after orthogonalisation against the (also constant) first one it is pure f32
/// rounding residue, which the power iteration then amplifies. Reproducing Fortran's
/// cuts therefore needs bit-exact f32 arithmetic in schmidt and the iteration.
fn gfortran_ran(flag: i32) -> f32 {
    const A: u64 = 16807;
    const M: u64 = 2_147_483_647;
    let seed: u64 = match flag {
        0 => unimplemented!("RAN(0) continues the libgfortran sequence; puu never passes 0"),
        1 => 123_459_876, // srand(0)
        f => f as u32 as u64,
    };
    let s = (A * seed) % M;
    (((s - 1) as u32) & !0x1FF) as f32 * (1.0 / 2_147_483_648.0)
}

fn schmidt(index: usize, x: &mut [Vec<f32>], xc: &mut [f32]) -> f32 {
    let n = x[index].len();
    let mut dotp = vec![0.0f32; index];
    for i in 0..n {
        for j in 0..index {
            dotp[j] += x[index][i] * x[j][i];
        }
    }
    let mut ssq = 1.0e-9f32;
    for i in 0..n {
        let mut e = x[index][i];
        for j in 0..index {
            e -= dotp[j] * x[j][i];
        }
        xc[i] = e;
        ssq += e * e;
    }
    let norm = ssq.sqrt();
    for i in 0..n {
        x[index][i] = xc[i] / norm;
    }
    norm
}

/// Second eigenvector of the contact correspondence matrix (index 1 of 0..2).
fn correspondence(n: usize, point: &[[usize; 2]], occx: &[f32], seed0: i32) -> Vec<f32> {
    let eigenvalues = 2;
    let mut x = vec![vec![0.0f32; n]; eigenvalues];
    let mut y = vec![vec![0.0f32; n]; eigenvalues];
    let mut rsums = vec![0.0f32; n];
    let mut csums = vec![0.0f32; n];
    let mut xc = vec![0.0f32; n];
    let mut f = vec![0.0f32; n];
    for (p, &o) in point.iter().zip(occx) {
        rsums[p[1] - 1] += o;
        csums[p[0] - 1] += o;
    }
    let mut prev_norm = 0.0f32;
    let mut norm;
    for index in 0..eigenvalues {
        for i in 0..n {
            x[index][i] = if index == 0 { 0.5 } else { gfortran_ran(seed0) };
        }
        norm = schmidt(index, &mut x, &mut xc);
        let mut round = 0;
        while (norm - prev_norm).abs() > 0.00001f32 * norm && round < 100 {
            round += 1;
            f.iter_mut().for_each(|v| *v = 0.0);
            for (p, &o) in point.iter().zip(occx) {
                f[p[1] - 1] += x[index][p[0] - 1] * o;
            }
            for j in 0..n {
                y[index][j] = f[j] / csums[j];
            }
            f.iter_mut().for_each(|v| *v = 0.0);
            for (p, &o) in point.iter().zip(occx) {
                f[p[0] - 1] += y[index][p[1] - 1] * o;
            }
            for i in 0..n {
                x[index][i] = f[i] / rsums[i];
            }
            prev_norm = norm;
            norm = schmidt(index, &mut x, &mut xc);
        }
    }
    x.swap_remove(1)
}

/// definesets: cuts where the eigenvector crosses `zerolevel`; returns icut[0..k]
/// (1-based residue positions).
fn definesets(x: &[f32], minseglen: usize, ltest: &[bool], zerolevel: f32, choice: u8) -> Vec<usize> {
    let nres = x.len();
    let set: Vec<u8> = x
        .iter()
        .map(|&v| if (choice == 1 && v < zerolevel) || (choice == 2 && v > zerolevel) { 2 } else { 1 })
        .collect();
    let mut icut: Vec<usize> = Vec::new();
    for i in 1..nres {
        if set[i] != set[i - 1] {
            icut.push(i); // Fortran icut(n)=i-1 with 1-based i
        }
    }
    if icut.last() == Some(&nres) {
        icut.pop();
    }
    let n = icut.len();
    for i in 0..n {
        let mut l = icut[i] as i64;
        if i > 0 {
            let mut j = i as i64 - 1;
            while j >= 0 && icut[j as usize] == 0 {
                j -= 1;
            }
            if j >= 0 {
                l = icut[i] as i64 - icut[j as usize] as i64;
            }
        }
        // ltest(icut(i)): icut(i) may already be 0 here (zeroed by i+1? no: only i-1)
        let ok = icut[i] >= 1 && ltest[icut[i] - 1];
        if l < minseglen as i64 || !ok {
            icut[i] = 0;
            if i > 0 {
                icut[i - 1] = 0;
            }
        }
    }
    let mut k: Vec<usize> = icut.into_iter().filter(|&c| c > 0).collect();
    if k.len() == 1 && k[0] == nres {
        k.clear();
    }
    k
}

/// tau-squared for cuts icut[0..bottom] (icut[bottom-1] = nres), using the
/// fastcont lower-triangle cumulative contact matrix.
fn tau(icut: &[usize], cont: &[f32], nres: usize, ncum: &[i32]) -> f32 {
    let bottom = icut.len();
    let c = |i: usize, j: usize| cont[(i - 1) * nres + (j - 1)];
    let mut x = 0.0f32;
    for iseg in 0..bottom {
        let mut jseg = iseg + 1;
        while jseg < bottom {
            for i in icut[jseg - 1] + 1..=icut[jseg] {
                x += c(i, icut[iseg]);
                if iseg > 0 {
                    x -= c(i, icut[iseg - 1]);
                }
            }
            jseg += 2;
        }
    }
    finish_tau(x, icut, ncum, nres)
}

fn tauslow(icut: &[usize], point: &[[usize; 2]], occx: &[f32], nres: usize, ncum: &[i32]) -> f32 {
    let bottom = icut.len();
    let mut set = vec![1u8; nres + 1];
    let mut iseg = 1;
    while iseg < bottom {
        for r in icut[iseg - 1] + 1..=icut[iseg] {
            set[r] = 2;
        }
        iseg += 2;
    }
    let mut x = 0.0f32;
    for (p, &o) in point.iter().zip(occx) {
        if set[p[0]] != set[p[1]] {
            x += o;
        }
    }
    finish_tau(x, icut, ncum, nres)
}

fn finish_tau(x: f32, icut: &[usize], ncum: &[i32], nres: usize) -> f32 {
    let x = x.max(1e-9);
    let mut m = 0.0f32;
    let mut jseg = 1;
    while jseg < icut.len() {
        m += (ncum[icut[jseg]] - ncum[icut[jseg - 1]]) as f32;
        jseg += 2;
    }
    let ntot = ncum[nres];
    let taux = m * (ntot as f32 - m) / (ntot.max(1) as f32) / x;
    taux.min(999.999)
}

/// fastcont: raw contact counts, then row-cumulative up to the diagonal.
fn fastcont(nres: usize, point: &[[usize; 2]], occx: &[f32]) -> Vec<f32> {
    let mut cont = vec![0.0f32; nres * nres];
    for (p, &o) in point.iter().zip(occx) {
        cont[(p[0] - 1) * nres + (p[1] - 1)] += o;
    }
    for i in 0..nres {
        for j in 1..=i {
            cont[i * nres + j] = cont[i * nres + j - 1] + cont[i * nres + j];
        }
    }
    cont
}

/// checkchain: residues within minseglen-1 of a chain end or break may not be cuts.
fn checkchain(ca: &[[f32; 3]], breakdist: f32, minseglen: usize) -> Vec<bool> {
    let nres = ca.len();
    let mut lbreak = vec![false; nres + 1]; // index 0..nres
    lbreak[0] = true;
    for i in 1..nres {
        lbreak[i] = distance(&ca[i - 1], &ca[i]) > breakdist;
    }
    lbreak[nres] = true;
    let mut ltest = vec![true; nres];
    let m = minseglen as i64;
    for i in 0..=nres {
        if lbreak[i] {
            for j in (i as i64 - m + 1)..=(i as i64 + m - 1) {
                if j >= 1 && j <= nres as i64 {
                    ltest[(j - 1) as usize] = false;
                }
            }
            if i >= 1 {
                ltest[i - 1] = true;
            }
        }
    }
    ltest
}

struct Work {
    nodes: Vec<PuuNode>, // 1-based via index-1
    set: Vec<usize>,     // per residue (1-based), node id
    stack: Vec<usize>,
}

fn new_node() -> PuuNode {
    PuuNode { child: [0, 0], node_type: ' ', tau: 0.0, tert: 0.0, nbp: 0, size: 0, ranges: Vec::new() }
}

/// dostack: binary decomposition until no set can be cut.
fn dostack(
    w: &mut Work,
    ca0: &[[f32; 3]],
    natom0: &[i32],
    contacts: &Contacts,
    p: &PuuParams,
) {
    let nres0 = ca0.len();
    while let Some(iset) = w.stack.pop() {
        // selectset
        let mut resix = Vec::new();
        let mut serix = vec![0usize; nres0 + 1];
        for i in 1..=nres0 {
            if w.set[i] == iset {
                resix.push(i);
                serix[i] = resix.len();
            }
        }
        let nres = resix.len();
        let ca: Vec<[f32; 3]> = resix.iter().map(|&i| ca0[i - 1]).collect();
        let mut point = Vec::new();
        let mut occx = Vec::new();
        for (pt, &o) in contacts.point.iter().zip(&contacts.occx) {
            if w.set[pt[0]] == iset && w.set[pt[1]] == iset {
                point.push([serix[pt[0]], serix[pt[1]]]);
                occx.push(o);
            }
        }
        let cont = if nres <= 500 { fastcont(nres, &point, &occx) } else { Vec::new() };
        let mut ncum = vec![0i32; nres + 1];
        for ix in 1..=nres {
            ncum[ix] = ncum[ix - 1] + natom0[resix[ix - 1]];
        }
        let ltest = checkchain(&ca, p.breakdist, p.minseglen);

        // eigen
        let x = correspondence(nres, &point, &occx, p.seed);
        let mut taubest = 0.0f32;
        let mut bestcut: Vec<usize> = Vec::new();
        for ires in 0..nres {
            if !ltest[ires] {
                continue;
            }
            for choice in [1u8, 2u8] {
                let mut icut = definesets(&x, p.minseglen, &ltest, x[ires], choice);
                if icut.is_empty() {
                    continue;
                }
                icut.push(nres);
                let t = if nres <= 500 {
                    tau(&icut, &cont, nres, &ncum)
                } else {
                    tauslow(&icut, &point, &occx, nres, &ncum)
                };
                if t > taubest {
                    taubest = t;
                    icut.pop();
                    bestcut = icut;
                }
            }
        }

        if taubest > 0.0 {
            savecut(w, iset, &bestcut, nres, &resix, taubest, p.minseglen);
        }
    }
}

fn savecut(w: &mut Work, iset: usize, bestcut: &[usize], nres: usize, resix: &[usize], taubest: f32, minseglen: usize) {
    let ndom = w.nodes.len();
    let bottom = bestcut.len();
    let mut cut: Vec<usize> = bestcut.to_vec(); // cut[k-1] = bestcut(k)
    cut.push(nres); // bestcut(bottom+1)=nres
    let bc = |k: usize| cut[k - 1]; // 1-based accessor

    w.nodes[iset - 1].tau = taubest;
    w.nodes[iset - 1].child = [ndom + 1, ndom + 2];
    w.nodes.push(new_node());
    w.nodes.push(new_node());

    // child 2: even segments
    let mut n2 = 0;
    let mut ranges2 = Vec::new();
    let mut iseg = 2;
    while iseg <= bottom + 1 {
        let j = bc(iseg - 1) + 1;
        for i in j..=bc(iseg) {
            w.set[resix[i - 1]] = ndom + 2;
            n2 += 1;
        }
        ranges2.push((resix[j - 1], resix[bc(iseg) - 1]));
        iseg += 2;
    }
    w.nodes[ndom + 1].ranges = ranges2;
    if n2 >= 2 * minseglen.max(1) {
        w.stack.push(ndom + 2);
    }

    // child 1: odd segments
    let mut n1 = 0;
    let mut ranges1 = Vec::new();
    let mut iseg = 1;
    while iseg <= bottom + 1 {
        let j = if iseg == 1 { 1 } else { bc(iseg - 1) + 1 };
        for i in j..=bc(iseg) {
            w.set[resix[i - 1]] = ndom + 1;
            n1 += 1;
        }
        ranges1.push((resix[j - 1], resix[bc(iseg) - 1]));
        iseg += 2;
    }
    w.nodes[ndom].ranges = ranges1;
    if n1 >= 2 * minseglen.max(1) {
        w.stack.push(ndom + 1);
    }
}

// ---------------------------------------------------------------------------
// Tree post-processing: patch, contstren, betasheets, findstrucdoms, pret2
// ---------------------------------------------------------------------------

/// patch: node ranges become the residues not claimed by a later (deeper) node,
/// as contiguous runs (at most MAXNSEG; Fortran overwrites the last one beyond).
fn patch(nodes: &mut [PuuNode], nres: usize) -> Vec<usize> {
    let ndom = nodes.len();
    let mut parent = vec![0usize; ndom + 1];
    for i in 1..=ndom {
        for &c in &nodes[i - 1].child {
            if c > 0 {
                parent[c] = i;
            }
        }
    }
    let mut resix = vec![1usize; nres + 1];
    for i in 2..=ndom {
        let p = parent[i];
        for &(a, b) in &nodes[i - 1].ranges.clone() {
            for k in a..=b {
                if resix[k] == p {
                    resix[k] = i;
                }
            }
        }
        let mut r: Vec<(usize, usize)> = Vec::new();
        let mut np = 0usize;
        let set_start = |r: &mut Vec<(usize, usize)>, np: &mut usize, s: usize| {
            *np += 1;
            if *np > MAXNSEG {
                *np = MAXNSEG;
            }
            if r.len() < *np {
                r.push((s, 0));
            } else {
                r[*np - 1].0 = s;
            }
        };
        if resix[1] == i {
            set_start(&mut r, &mut np, 1);
        }
        for j in 1..nres {
            if resix[j] != i && resix[j + 1] == i {
                set_start(&mut r, &mut np, j + 1);
            }
            if resix[j] == i && resix[j + 1] != i {
                r[np - 1].1 = j;
            }
        }
        if resix[nres] == i {
            r[np - 1].1 = nres;
        }
        r.truncate(np);
        nodes[i - 1].ranges = r;
    }
    nodes[0].size = nres;
    for i in 2..=ndom {
        nodes[i - 1].size = nodes[i - 1].ranges.iter().map(|&(a, b)| b + 1 - a).sum();
    }
    parent
}

/// contstren with natom=1 per domain (see module docs).
fn contstren(nodes: &mut [PuuNode], contacts: &Contacts, nres: usize) {
    const FOUR: usize = 4;
    if nres <= 500 {
        let cont = fastcont(nres, &contacts.point, &contacts.occx);
        let c = |i: usize, j: usize| cont[(i - 1) * nres + (j - 1)];
        for node in nodes.iter_mut() {
            let mut x = 1e-9f32;
            for &(k, l) in &node.ranges {
                let mut j = k + FOUR;
                while j <= l {
                    x += c(j, j - FOUR);
                    if k > 1 {
                        x -= c(j, k - 1);
                    }
                    j += 1;
                }
            }
            node.tert = x / 1.0 * 100.0;
        }
    } else {
        for node in nodes.iter_mut() {
            let mut set = vec![2u8; nres + 1];
            for &(a, b) in &node.ranges {
                for r in a..=b {
                    set[r] = 1;
                }
            }
            let mut x = 1e-9f32;
            for (p, &o) in contacts.point.iter().zip(&contacts.occx) {
                if (p[0] as i64 - p[1] as i64).abs() >= FOUR as i64 && set[p[0]] == 1 && set[p[1]] == 1 {
                    x += o;
                }
            }
            node.tert = x / 1.0 * 100.0;
        }
    }
}

fn betasheets(nodes: &mut [PuuNode], bp: &[[i32; 2]], nres: usize) {
    let snapshot: Vec<PuuNode> = nodes.to_vec();
    for node in nodes.iter_mut() {
        let [j, k] = node.child;
        node.nbp = if j > 0 && k > 0 {
            let mut idom = vec![false; nres + 1];
            for &(a, b) in &snapshot[j - 1].ranges {
                for r in a..=b {
                    idom[r] = true;
                }
            }
            let mut n = 0;
            for &(a, b) in &snapshot[k - 1].ranges {
                for r in a..=b {
                    let (k1, k2) = (bp[r][0], bp[r][1]);
                    if k1 > 0 && k2 > 0 {
                        if (k1 as usize) <= nres && idom[k1 as usize] {
                            n += 1;
                        }
                        if (k2 as usize) <= nres && idom[k2 as usize] {
                            n += 1;
                        }
                    }
                }
            }
            n
        } else {
            0
        };
    }
}

/// testsplit; returns lstop.
fn testsplit(nodes: &[PuuNode], i: usize, j: usize, k: usize, p: &PuuParams, minsize: usize) -> bool {
    let n = &nodes[i - 1];
    if j == 0 || k == 0 {
        true
    } else if n.size < minsize {
        true
    } else if n.tau > p.taucutoff {
        false
    } else if n.nbp > 0 && !p.lsplit {
        true
    } else if n.tert < p.compact {
        true
    } else if nodes[j - 1].tert <= p.compact && nodes[k - 1].tert <= p.compact {
        true
    } else if nodes[j - 1].tert > p.compact && nodes[k - 1].tert > p.compact {
        false
    } else {
        true
    }
}

fn delete_subtree(nodes: &mut [PuuNode], del: &mut [bool], knode: usize) {
    let mut stack = Vec::new();
    if knode > 0 {
        stack.push(knode);
    }
    while let Some(inode) = stack.pop() {
        nodes[inode - 1].node_type = '-';
        del[inode] = true;
        for &c in &nodes[inode - 1].child.clone() {
            if c > 0 {
                stack.push(c);
            }
        }
    }
}

fn findstrucdoms(nodes: &mut [PuuNode], p: &PuuParams, minsize: usize) {
    let ndom = nodes.len();
    let halfminsize = minsize / 2;
    let mut del = vec![false; ndom + 1];
    for n in nodes.iter_mut() {
        n.node_type = '*';
    }
    for i in 1..=ndom {
        let [mut j, k] = nodes[i - 1].child;
        if del[i] {
            j = 0;
        }
        let mut lstop = testsplit(nodes, i, j, k, p, minsize);
        // 1nrd hack
        if j > 0 && k > 0 && lstop {
            let (m, n) = if nodes[k - 1].size < nodes[j - 1].size { (k, j) } else { (j, k) };
            let mut lstop1 = true;
            if nodes[m - 1].size < halfminsize
                && nodes[n - 1].size > minsize
                && nodes[m - 1].tert < p.compact
                && nodes[n - 1].tert > p.compact
            {
                let [c1, c2] = nodes[n - 1].child;
                lstop1 = testsplit(nodes, n, c1, c2, p, minsize);
            }
            if !lstop1 {
                lstop = false;
            }
        }
        if lstop {
            delete_subtree(nodes, &mut del, j);
            delete_subtree(nodes, &mut del, k);
        } else {
            del[i] = true;
            nodes[i - 1].node_type = '+';
        }
    }
}

/// Run puu on one chain: the tree puu writes to units.puu.
pub fn puu(s: &PuuStructure, dssp: &DsspInput, p: &PuuParams) -> Vec<PuuNode> {
    let nres = s.nres;

    // grepatom 'CA ' (last one per residue wins), on uncentered coordinates
    let mut ca = vec![[0.0f32; 3]; nres];
    for a in &s.atoms {
        if a.res >= 1 && a.res <= nres && &a.name == b"CA " {
            ca[a.res - 1] = a.xyz;
        }
    }
    let mut natom0 = vec![0i32; nres + 1];
    for a in &s.atoms {
        if a.res >= 1 && a.res <= nres {
            natom0[a.res] += 1;
        }
    }

    let mut contacts = get_contacts(s);
    let bp = add_hbonds(&mut contacts, dssp, nres, p.hbond);

    // initialize_one + dostack
    let mut root = new_node();
    root.ranges = vec![(1.min(nres), nres)];
    let mut w = Work { nodes: vec![root], set: vec![0; nres + 1], stack: vec![1] };
    for i in 1..=nres {
        w.set[i] = 1;
    }
    dostack(&mut w, &ca, &natom0, &contacts, p);
    let mut nodes = w.nodes;

    let parent = patch(&mut nodes, nres);
    contstren(&mut nodes, &contacts, nres);
    betasheets(&mut nodes, &bp, nres);
    findstrucdoms(&mut nodes, p, 80);

    // pret2: structural domains of <= 40 residues become '='
    if nodes[0].size == 0 {
        return Vec::new();
    }
    for n in nodes.iter_mut() {
        if (n.node_type == '*' || n.node_type == '=') && n.size <= 40 {
            n.node_type = '=';
        }
    }

    // writetree1: keep nodes down to printdepth, renumber
    let ndom = nodes.len();
    let mut dep = vec![0i32; ndom + 1];
    let mut keep = Vec::new();
    let mut tsil = vec![0usize; ndom + 1];
    for i in 1..=ndom {
        dep[i] = match nodes[i - 1].node_type {
            '+' => 1,
            '-' => {
                let pi = parent[i];
                let mut d = dep[pi] - 1;
                if pi > 0 && nodes[pi - 1].node_type == '=' {
                    d -= 1;
                }
                d
            }
            _ => 0,
        };
        if dep[i] >= p.printdepth {
            keep.push(i);
            tsil[i] = keep.len();
        }
    }
    keep.iter()
        .map(|&j| {
            let mut n = nodes[j - 1].clone();
            n.child = if n.child[0] > 0 { [tsil[n.child[0]], tsil[n.child[1]]] } else { [0, 0] };
            n
        })
        .collect()
}

/// The units.puu tree as residue-level DomainNodes (input to `io::puutos`).
pub fn to_domain_tree(nodes: &[PuuNode]) -> Vec<crate::types::DomainNode> {
    nodes.iter().enumerate()
        .map(|(i, n)| crate::types::DomainNode {
            index: i + 1,
            node_type: crate::types::NodeType::from_char(n.node_type),
            left_child: n.child[0],
            right_child: n.child[1],
            nseg: n.ranges.len(),
            segments: n.ranges.iter().map(|&(a, b)| (a as u32, b as u32)).collect(),
        })
        .collect()
}

/// Parse the NUNITS block of a units.puu file (for validation against Fortran).
pub fn read_units_puu<P: AsRef<Path>>(path: P) -> Result<Vec<PuuNode>, PuuError> {
    let text = read_text(path)?;
    let mut lines = text.lines().skip_while(|l| !l.starts_with("NUNITS"));
    let header = lines.next().ok_or_else(|| PuuError::Io("no NUNITS".into()))?;
    let n: usize = header[10..].trim().parse().map_err(|_| PuuError::Io("bad NUNITS".into()))?;
    lines.next(); // column header
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        // (i5,1x,2i4,1x,a1,a32,2f9.3,i8,2i7)
        let l = fixed(lines.next().ok_or_else(|| PuuError::Io("truncated".into()))?, 90);
        let s = |a: usize, b: usize| String::from_utf8_lossy(&l[a..b]).trim().to_string();
        let c1: usize = s(6, 10).parse().unwrap_or(0);
        let c2: usize = s(10, 14).parse().unwrap_or(0);
        let node_type = l[15] as char;
        let tau: f32 = s(48, 57).parse().unwrap_or(0.0);
        let tert: f32 = s(57, 66).parse().unwrap_or(0.0);
        let nbp: i32 = s(66, 74).parse().unwrap_or(0);
        let size: usize = s(74, 81).parse().unwrap_or(0);
        let nseg: usize = s(81, 88).parse().unwrap_or(0);
        let mut ranges = Vec::with_capacity(nseg);
        for _ in 0..nseg {
            let r = lines.next().ok_or_else(|| PuuError::Io("truncated".into()))?;
            let v: Vec<usize> = r.split_whitespace().take(2).filter_map(|x| x.parse().ok()).collect();
            ranges.push((v[0], v[1]));
        }
        out.push(PuuNode { child: [c1, c2], node_type, tau, tert, nbp, size, ranges });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::gfortran_ran;

    /// Bit patterns printed by gfortran -O3 for `x=ran(s)` (see gfortran_ran).
    #[test]
    fn test_gfortran_ran_matches_libgfortran() {
        assert_eq!(gfortran_ran(1_234_567_811).to_bits(), 0x3E3BF860);
        assert_eq!(gfortran_ran(7).to_bits(), 0x38650000);
    }
}
