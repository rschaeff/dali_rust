//! DSSP secondary structure — port of DaliLite's dsspcmbi (Kabsch & Sander 1983,
//! CMBI version by E. Krieger, p2c translation in src/dsspcmbi/DsspCMBI.c).
//!
//! Ports what DALI's import consumes: the residue list with chain-break ('!')
//! entries, H-bonds, bridges, ladders, sheets, bridge partners (BP1/BP2), turns,
//! helices, bends and the summary structure column. Accessibility and the
//! disulfide relabelling of cysteines are not ported (nothing downstream reads
//! them). Arithmetic is f64 in the C evaluation order; the H-bond energy is
//! rounded with floor(x+0.5) exactly as the C does.
//!
//! Residue indices follow the C `chain[]` array: 1..=lchain, break entries
//! included, index 0 unused.

use std::path::Path;

const DIST: f64 = 0.5;
const BREAKDIST: f64 = 2.5;
const CADIST: f64 = 9.0;
const Q: f64 = -27888.0;
const HBLOW: i64 = -9900;
const HBHIGH: i64 = -500;
const RADIAN: f64 = 57.29578;
const PIHALF: f64 = 1.570796;
const PI: f64 = 3.141593;
const TWOPI: f64 = 6.283185;
const MAXBRIDGE: usize = 1000;

// ss[] column indices (the C `structure` enum)
const SYMBOL: usize = 0;
const TURN3: usize = 1;
const BEND: usize = 4;
const CHIRALITY: usize = 5;
const BETA1: usize = 6;

type V3 = [f64; 3];

#[derive(Debug, Clone, Copy, Default)]
struct HBond {
    residue: usize,
    energy: i64,
}

/// One entry of the C `chain[]` array.
#[derive(Debug, Clone)]
pub struct DsspResidue {
    /// resnum (4) + insertion code + chain, as in DSSP columns 7-12
    pub aaident: [u8; 6],
    /// one-letter code; b'!' for a chain-break entry
    pub aa: u8,
    pub threeletter: [u8; 4],
    /// ss[0]=summary, 1..3 turn3..5, 4 bend, 5 chirality, 6..7 beta1/beta2 ladder
    pub ss: [u8; 8],
    /// bridge partners (DSSP numbering, 0 = none)
    pub partner: [i64; 2],
    pub sheetlabel: u8,
    pub kappa: f64,
    pub alpha: f64,
    acceptor: [HBond; 2],
    donor: [HBond; 2],
    pub n: V3,
    pub ca: V3,
    pub c: V3,
    pub o: V3,
    pub h: V3,
}

impl DsspResidue {
    fn zero() -> Self {
        DsspResidue {
            aaident: [b' '; 6],
            aa: b'!',
            threeletter: *b"    ",
            ss: [b' '; 8],
            partner: [0, 0],
            sheetlabel: b' ',
            kappa: 360.0,
            alpha: 360.0,
            acceptor: [HBond::default(); 2],
            donor: [HBond::default(); 2],
            n: [0.0; 3],
            ca: [0.0; 3],
            c: [0.0; 3],
            o: [0.0; 3],
            h: [0.0; 3],
        }
    }

    pub fn is_break(&self) -> bool {
        self.aa == b'!'
    }

    pub fn chain_id(&self) -> u8 {
        self.aaident[5]
    }
}

/// dsspcmbi result: `res[1..=lchain]`.
#[derive(Debug, Clone)]
pub struct Dssp {
    pub res: Vec<DsspResidue>,
}

#[derive(Debug)]
pub enum DsspError {
    Io(String),
    NoResidues,
}

impl std::fmt::Display for DsspError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DsspError::Io(s) => write!(f, "{}", s),
            DsspError::NoResidues => write!(f, "no residue with complete backbone"),
        }
    }
}

impl std::error::Error for DsspError {}

#[inline]
fn distsq(u: &V3, v: &V3) -> f64 {
    (u[0] - v[0]) * (u[0] - v[0]) + (u[1] - v[1]) * (u[1] - v[1]) + (u[2] - v[2]) * (u[2] - v[2])
}

#[inline]
fn distance(u: &V3, v: &V3) -> f64 {
    distsq(u, v).sqrt()
}

#[inline]
fn diff(x: &V3, y: &V3) -> V3 {
    [x[0] - y[0], x[1] - y[1], x[2] - y[2]]
}

#[inline]
fn dot(x: &V3, y: &V3) -> f64 {
    x[0] * y[0] + x[1] * y[1] + x[2] * y[2]
}

#[inline]
fn cross(x: &V3, y: &V3) -> V3 {
    [
        x[1] * y[2] - y[1] * x[2],
        x[2] * y[0] - y[2] * x[0],
        x[0] * y[1] - y[0] * x[1],
    ]
}

fn atan2_c(y: f64, x: f64) -> f64 {
    let mut z = if x != 0.0 {
        (y / x).atan()
    } else if y > 0.0 {
        PIHALF
    } else if y < 0.0 {
        -PIHALF
    } else {
        TWOPI
    };
    if x >= 0.0 {
        return z;
    }
    if y > 0.0 {
        z += PI;
    } else {
        z -= PI;
    }
    z
}

fn dihedral(v1: &V3, v2: &V3, v3: &V3, v4: &V3) -> f64 {
    let v12 = diff(v1, v2);
    let v43 = diff(v4, v3);
    let z = diff(v2, v3);
    let p = cross(&z, &v12);
    let x = cross(&z, &v43);
    let y = cross(&z, &x);
    let mut u = dot(&x, &x);
    let mut v = dot(&y, &y);
    if u <= 0.0 || v <= 0.0 {
        return 360.0;
    }
    u = dot(&p, &x) / u.sqrt();
    v = dot(&p, &y) / v.sqrt();
    if u != 0.0 || v != 0.0 {
        return atan2_c(v, u) * RADIAN;
    }
    360.0
}

fn cosangle(v1: &V3, v2: &V3, v3: &V3, v4: &V3) -> f64 {
    let u = diff(v1, v2);
    let v = diff(v3, v4);
    let x = dot(&u, &u) * dot(&v, &v);
    if x > 0.0 { dot(&u, &v) / x.sqrt() } else { 0.0 }
}

// ---------------------------------------------------------------------------
// Input (Readcard / Putresidue / Getresidue / Inputcoordinates)
// ---------------------------------------------------------------------------

fn onelettercode(aaa: &[u8]) -> u8 {
    const AASYMBOL: &[u8; 50] = b"ARNDCEQGHILKMFPSTWYVBZXXXXXXXXXXXXXXXX--CCCCIPPPW-";
    const AMINOACID: &[u8; 150] = b"ALAARGASNASPCYSGLUGLNGLYHISILE\
LEULYSMETPHEPROSERTHRTRPTYRVAL\
ASXGLXACDALBALIABUAROBASBETHSE\
HYPHYLORNPCASARTAUTHYUNKACEFOR\
CYHCSHCSSCYXILUPRZPR0CPRTRYHOH";
    let mut k = 1;
    let mut i = 1;
    while k < 51 {
        if AMINOACID[i - 1] == aaa[0] && AMINOACID[i] == aaa[1] && AMINOACID[i + 1] == aaa[2] {
            let a = AASYMBOL[k - 1];
            return if a == b'-' { b'-' } else { a };
        }
        i += 3;
        k += 1;
    }
    b'X' // "Let's assume that anything is an aminoacid."
}

/// fscanf("%lf") x n from the stream position: strtod's longest valid prefix, so
/// run-together fields split where a number can no longer continue
/// ("640.6121103.3031323.682" -> 640.6121103, .3031323, .682 -- garbage, but it is
/// what dsspcmbi computes with for coordinates >= 1000 A).
fn scan_doubles(s: &[u8], n: usize) -> Vec<f64> {
    let mut out = Vec::with_capacity(n);
    let mut i = 0;
    while out.len() < n && i < s.len() {
        while i < s.len() && (s[i] as char).is_ascii_whitespace() {
            i += 1;
        }
        let start = i;
        if i < s.len() && (s[i] == b'-' || s[i] == b'+') {
            i += 1;
        }
        while i < s.len() && s[i].is_ascii_digit() {
            i += 1;
        }
        if i < s.len() && s[i] == b'.' {
            i += 1;
            while i < s.len() && s[i].is_ascii_digit() {
                i += 1;
            }
        }
        if i < s.len() && (s[i] == b'e' || s[i] == b'E') {
            let save = i;
            i += 1;
            if i < s.len() && (s[i] == b'-' || s[i] == b'+') {
                i += 1;
            }
            let d = i;
            while i < s.len() && s[i].is_ascii_digit() {
                i += 1;
            }
            if i == d {
                i = save;
            }
        }
        // Rust's parser rejects a bare "." or "-.5"-style forms strtod accepts
        let tok = std::str::from_utf8(&s[start..i]).unwrap_or("");
        let fixed = if tok.starts_with('.') { format!("0{}", tok) }
            else if tok.starts_with("-.") { format!("-0{}", &tok[1..]) }
            else if tok.starts_with("+.") { format!("0{}", &tok[1..]) }
            else { tok.to_string() };
        match fixed.parse::<f64>().ok() {
            Some(v) => out.push(v),
            None => break,
        }
    }
    out
}

struct Reader {
    chain: Vec<DsspResidue>, // index 0..=lchain (+1 slack)
    lchain: usize,
    resinfo: DsspResidue,
    nmissing: bool,
    camissing: bool,
    cmissing: bool,
    omissing: bool,
}

impl Reader {
    fn at(&self, i: usize) -> &DsspResidue {
        &self.chain[i]
    }

    fn ensure(&mut self, i: usize) {
        while self.chain.len() <= i + 1 {
            self.chain.push(DsspResidue::zero());
        }
    }

    fn nochainbreak(&self, i: usize, j: usize) -> bool {
        if !(i >= 1 && i <= j) {
            return false;
        }
        (i..=j).all(|k| k < self.chain.len() && !self.chain[k].is_break())
    }

    fn putresidue(&mut self) {
        let complete = !(self.nmissing || self.camissing || self.cmissing || self.omissing);
        let mut superposed = false;
        let r = &self.resinfo;
        if r.aaident[4] != b' ' && complete {
            for i in (1..=self.lchain).rev() {
                let c = self.at(i);
                if distance(&c.n, &r.n) < DIST || distance(&c.ca, &r.ca) < DIST || distance(&c.c, &r.c) < DIST {
                    superposed = true;
                    break;
                }
            }
        }
        if r.aaident[4] == b' ' && complete && self.lchain > 0 {
            let c = self.at(self.lchain);
            if distance(&c.n, &r.n) < DIST || distance(&c.ca, &r.ca) < DIST || distance(&c.c, &r.c) < DIST {
                superposed = true;
            }
        }
        if complete && !superposed {
            let mut res = self.resinfo.clone();
            res.h = res.n;
            if self.nochainbreak(self.lchain, self.lchain)
                && distance(&self.at(self.lchain).c, &res.n) > BREAKDIST
            {
                self.lchain += 1; // keep '!' at lchain
                self.ensure(self.lchain);
            }
            if self.nochainbreak(self.lchain, self.lchain) && res.aa != b'P' {
                let prev = self.at(self.lchain);
                let dco = distance(&prev.c, &prev.o);
                for i in 0..3 {
                    res.h[i] = res.n[i] + (prev.c[i] - prev.o[i]) / dco;
                }
            }
            self.lchain += 1;
            self.ensure(self.lchain);
            self.chain[self.lchain] = res;
        }
        if self.nochainbreak(self.lchain, self.lchain) && !complete {
            self.lchain += 1;
            self.ensure(self.lchain);
        }
        self.resinfo = DsspResidue::zero();
        self.nmissing = true;
        self.camissing = true;
        self.cmissing = true;
        self.omissing = true;
    }

    fn getresidue(&mut self, atomname: &[u8; 4], xyz: V3) {
        let hydrogen = (atomname[0].is_ascii_digit() || atomname[0] == b' ')
            && (atomname[1] == b'D' || atomname[1] == b'H');
        if hydrogen {
            return;
        }
        match atomname {
            b" N  " => {
                self.nmissing = false;
                self.resinfo.n = xyz;
            }
            b" CA " => {
                self.camissing = false;
                self.resinfo.ca = xyz;
            }
            b" C  " => {
                self.cmissing = false;
                self.resinfo.c = xyz;
            }
            b" O  " => {
                self.omissing = false;
                self.resinfo.o = xyz;
            }
            _ => {} // side chain: accessibility only
        }
    }
}

#[derive(PartialEq)]
enum Card {
    Atom,
    Ter,
    End,
    Other,
}

/// Inputcoordinates: PDB text -> chain[] with break entries.
fn input_coordinates(text: &str) -> Result<Reader, DsspError> {
    let mut rd = Reader {
        chain: vec![DsspResidue::zero(), DsspResidue::zero()],
        lchain: 0,
        resinfo: DsspResidue::zero(),
        nmissing: true,
        camissing: true,
        cmissing: true,
        omissing: true,
    };
    let mut lastchain: u8 = 0;
    let mut lastresno: i64 = -1;
    let mut altrename: u8 = 0;
    // Readcard's `key` keeps characters a short line does not overwrite
    let mut key = [0u8; 6];

    for raw in text.lines() {
        let b = raw.as_bytes();
        // skip to the first uppercase character
        let p = match b.iter().position(|c| c.is_ascii_uppercase()) {
            Some(p) => p,
            None => continue,
        };
        let line = &b[p..];
        for l in 0..6 {
            if l < line.len() {
                key[l] = line[l];
            }
        }
        let card = match &key {
            b"ATOM  " | b"HETATM" => Card::Atom,
            b"TER   " => Card::Ter,
            b"END   " | b"ENDMDL" => Card::End,
            b"MODEL " => {
                if lastchain != 0 { Card::End } else { Card::Other }
            }
            _ => Card::Other,
        };
        match card {
            Card::Atom => {
                let get = |i: usize| -> u8 {
                    match line.get(i) {
                        Some(&c) => c,
                        None => b' ',
                    }
                };
                let atomname = [get(12), get(13), get(14), get(15)];
                let mut altloc = get(16);
                let aaname = [get(17), get(18), get(19), b' '];
                let residuename = onelettercode(&aaname);
                let c = get(21);
                if c != lastchain {
                    lastchain = c;
                    lastresno = -1;
                    altrename = 0;
                }
                let resnum = [get(22), get(23), get(24), get(25)];
                let n: i64 = {
                    // atoi
                    let s = String::from_utf8_lossy(&resnum).to_string();
                    let t = s.trim_start();
                    let mut end = 0;
                    let tb = t.as_bytes();
                    if end < tb.len() && (tb[end] == b'-' || tb[end] == b'+') {
                        end += 1;
                    }
                    while end < tb.len() && tb[end].is_ascii_digit() {
                        end += 1;
                    }
                    t[..end].parse().unwrap_or(0)
                };
                let icode = get(26);
                let aaident = [resnum[0], resnum[1], resnum[2], resnum[3], icode, c];
                if n > lastresno {
                    lastresno = n;
                    altrename = if altloc > b'A' || (altloc > b'1' && altloc <= b'9') { altloc } else { 0 };
                }
                if altloc == altrename {
                    altloc = b'A';
                }
                let coords = if line.len() > 27 { scan_doubles(&line[27..], 3) } else { Vec::new() };
                let xyz = [
                    *coords.first().unwrap_or(&0.0),
                    *coords.get(1).unwrap_or(&0.0),
                    *coords.get(2).unwrap_or(&0.0),
                ];
                if residuename != b'-' && (altloc == b'A' || altloc == b'1' || altloc == b' ') {
                    if rd.resinfo.aaident != aaident || rd.resinfo.threeletter != aaname {
                        rd.putresidue();
                        rd.resinfo.aaident = aaident;
                        rd.resinfo.aa = residuename;
                        rd.resinfo.threeletter = aaname;
                    }
                    rd.getresidue(&atomname, xyz);
                }
            }
            Card::Ter => rd.putresidue(),
            Card::End => {
                rd.putresidue();
                break;
            }
            Card::Other => {}
        }
        if card == Card::End {
            break;
        }
    }
    // EOF without END: Readcard returns endcard
    if !(rd.resinfo.aa == b'!' && rd.nmissing && rd.camissing && rd.cmissing && rd.omissing) {
        rd.putresidue();
    }
    if !rd.nochainbreak(rd.lchain, rd.lchain) && rd.lchain > 0 {
        rd.lchain -= 1;
    }
    if rd.lchain < 1 {
        return Err(DsspError::NoResidues);
    }
    rd.chain.truncate(rd.lchain + 2); // keep one zero entry past the end
    Ok(rd)
}

// ---------------------------------------------------------------------------
// Secondary structure
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, PartialEq, Debug)]
enum BType {
    Parallel,
    Antiparallel,
}

#[derive(Clone, Debug)]
struct Bridge {
    sheetname: u8,
    laddername: u8,
    btyp: BType,
    linkset: Vec<bool>, // indexed 1..=nbridge
    ib: usize,
    ie: usize,
    jb: usize,
    je: usize,
    from: usize,
    towards: usize,
}

struct State {
    c: Vec<DsspResidue>,
    lchain: usize,
    bridges: Vec<Bridge>, // bridgetable[0..nbridge]
}

impl State {
    fn nochainbreak(&self, i: usize, j: usize) -> bool {
        if !(i >= 1 && i <= j) {
            return false;
        }
        (i..=j).all(|k| k < self.c.len() && !self.c[k].is_break())
    }

    fn testbond(&self, i: usize, j: usize) -> bool {
        let r = &self.c[i];
        (r.acceptor[0].residue == j && r.acceptor[0].energy < HBHIGH)
            || (r.acceptor[1].residue == j && r.acceptor[1].energy < HBHIGH)
    }

    fn flagchirality(&mut self) {
        let l = self.lchain as i64;
        for i in 2..=(l - 2).max(1) as usize {
            if i as i64 > l - 2 {
                break;
            }
            if self.nochainbreak(i - 1, i + 2) {
                let a = dihedral(&self.c[i - 1].ca, &self.c[i].ca, &self.c[i + 1].ca, &self.c[i + 2].ca);
                self.c[i].alpha = a;
                self.c[i].ss[CHIRALITY] = if a < 0.0 { b'-' } else { b'+' };
            }
        }
        let mut i = 3usize;
        while (i as i64) <= l - 2 {
            if self.nochainbreak(i - 2, i + 2) {
                let ckap = cosangle(&self.c[i].ca, &self.c[i - 2].ca, &self.c[i + 2].ca, &self.c[i].ca);
                let skap = (1.0 - ckap * ckap).sqrt();
                self.c[i].kappa = RADIAN * atan2_c(skap, ckap);
            }
            i += 1;
        }
    }

    fn bondenergy(&self, i: usize, j: usize) -> i64 {
        let d = &self.c[i];
        if d.aa == b'P' {
            return 0;
        }
        let a = &self.c[j];
        let dho = distance(&d.h, &a.o);
        let dhc = distance(&d.h, &a.c);
        let dnc = distance(&d.n, &a.c);
        let dno = distance(&d.n, &a.o);
        let hbe = if dho < DIST || dhc < DIST || dnc < DIST || dno < DIST {
            HBLOW
        } else {
            (Q / dho - Q / dhc + Q / dnc - Q / dno + 0.5).floor() as i64
        };
        if hbe > HBLOW { hbe } else { HBLOW }
    }

    fn setbonds(&mut self, i: usize, j: usize) {
        let e = self.bondenergy(i, j);
        let upd = |b: &mut [HBond; 2], hb: HBond| {
            if hb.energy < b[0].energy {
                b[1] = b[0];
                b[0] = hb;
            } else if hb.energy < b[1].energy {
                b[1] = hb;
            }
        };
        upd(&mut self.c[i].acceptor, HBond { residue: j, energy: e });
        upd(&mut self.c[j].donor, HBond { residue: i, energy: e });
    }

    fn flaghydrogenbonds(&mut self) {
        for i in 1..=self.lchain {
            if !self.nochainbreak(i, i) {
                continue;
            }
            for j in i + 1..=self.lchain {
                if self.nochainbreak(j, j) && distance(&self.c[i].ca, &self.c[j].ca) < CADIST {
                    self.setbonds(i, j);
                    if j != i + 1 {
                        self.setbonds(j, i);
                    }
                }
            }
        }
    }

    fn ladder(&mut self, i: usize, j: usize, b: BType) {
        if i >= j {
            return;
        }
        for k in 0..self.bridges.len() {
            let w = &self.bridges[k];
            let found = w.btyp == b
                && i == w.ie + 1
                && self.nochainbreak(w.ie, i)
                && ((j == w.je + 1 && b == BType::Parallel && self.nochainbreak(w.je, j))
                    || (j + 1 == w.jb && b == BType::Antiparallel && self.nochainbreak(j, w.jb)));
            if found {
                let w = &mut self.bridges[k];
                w.ie += 1;
                if b == BType::Parallel {
                    w.je += 1;
                } else {
                    w.jb -= 1;
                }
                return;
            }
        }
        assert!(self.bridges.len() < MAXBRIDGE, "bridgetable overflow");
        self.bridges.push(Bridge {
            sheetname: 0,
            laddername: 0,
            btyp: b,
            linkset: Vec::new(),
            ib: i,
            ie: i,
            jb: j,
            je: j,
            from: 0,
            towards: 0,
        });
    }

    fn testbridge(&mut self, i: usize) {
        if !self.nochainbreak(i - 1, i + 1) {
            return;
        }
        let (mut j1, mut j2) = (0usize, 0usize);
        let mut j = i + 3;
        while j2 == 0 && j < self.lchain {
            if self.nochainbreak(j - 1, j + 1) {
                let b = if (self.testbond(i + 1, j) && self.testbond(j, i - 1))
                    || (self.testbond(j + 1, i) && self.testbond(i, j - 1))
                {
                    Some(BType::Parallel)
                } else if (self.testbond(i + 1, j - 1) && self.testbond(j + 1, i - 1))
                    || (self.testbond(j, i) && self.testbond(i, j))
                {
                    Some(BType::Antiparallel)
                } else {
                    None
                };
                if let Some(b) = b {
                    if j1 == 0 {
                        j1 = j;
                        self.ladder(i, j, b);
                    } else if j != j1 {
                        j2 = j;
                        self.ladder(i, j, b);
                    }
                }
            }
            j += 1;
        }
    }

    fn extendladder(&mut self) {
        let nb = self.bridges.len();
        for i in 1..=nb {
            let mut j = i + 1;
            while j <= nb && self.bridges[i - 1].towards == 0 {
                let (wib, wie, wjb, wje, wbt) = {
                    let w = &self.bridges[i - 1];
                    (w.ib, w.ie, w.jb, w.je, w.btyp)
                };
                let _ = wib;
                let (ib1, jb1, je1, bt1, from1) = {
                    let o = &self.bridges[j - 1];
                    (o.ib, o.jb, o.je, o.btyp, o.from)
                };
                let d = |a: usize, b: usize| a as i64 - b as i64;
                let mut bulge = self.nochainbreak(wie, ib1) && d(ib1, wie) < 6 && bt1 == wbt && from1 == 0;
                if bulge {
                    bulge = match wbt {
                        BType::Parallel => {
                            ((d(jb1, wje) < 6 && d(ib1, wie) < 3) || d(jb1, wje) < 3) && self.nochainbreak(wje, jb1)
                        }
                        BType::Antiparallel => {
                            ((d(wjb, je1) < 6 && d(ib1, wie) < 3) || d(wjb, je1) < 3) && self.nochainbreak(je1, wjb)
                        }
                    };
                }
                if bulge {
                    self.bridges[i - 1].towards = j;
                    self.bridges[j - 1].from = i;
                }
                j += 1;
            }
        }
        for i in 1..=nb {
            if self.bridges[i - 1].from == 0 {
                let mut set = vec![false; nb + 1];
                let mut j = i;
                loop {
                    set[j] = true;
                    j = self.bridges[j - 1].towards;
                    if j == 0 {
                        break;
                    }
                }
                self.bridges[i - 1].linkset = set.clone();
                let mut j = self.bridges[i - 1].towards;
                while j != 0 {
                    self.bridges[j - 1].linkset = set.clone();
                    j = self.bridges[j - 1].towards;
                }
            }
        }
    }

    fn link(&self, l1: usize, l2: usize) -> bool {
        let a = &self.bridges[l1 - 1];
        let b = &self.bridges[l2 - 1];
        (a.ie >= b.ib && a.ib <= b.ie)
            || (a.ie >= b.jb && a.ib <= b.je)
            || (a.je >= b.ib && a.jb <= b.ie)
            || (a.je >= b.jb && a.jb <= b.je)
    }

    fn findsheet(&self, ladderset: &mut Vec<bool>) -> Vec<bool> {
        let nb = self.bridges.len();
        let mut sheetset = vec![false; nb + 1];
        let l1 = match (1..=nb).find(|&l| ladderset[l]) {
            Some(l) => l,
            None => return sheetset,
        };
        sheetset = self.bridges[l1 - 1].linkset.clone();
        loop {
            let mut finish = true;
            for l1 in 1..=nb {
                if !sheetset[l1] {
                    continue;
                }
                for l2 in 1..=nb {
                    if ladderset[l2] && self.link(l1, l2) {
                        let ls = &self.bridges[l2 - 1].linkset;
                        for k in 1..=nb {
                            if ls[k] {
                                sheetset[k] = true;
                                ladderset[k] = false;
                            }
                        }
                        finish = false;
                    }
                }
            }
            if finish {
                break;
            }
        }
        sheetset
    }

    fn sheet(&mut self) {
        let nb = self.bridges.len();
        let mut ladderset = vec![true; nb + 1];
        ladderset[0] = false;
        let mut ccs: u8 = b'A' - 1;
        let mut asci: u8 = 64;
        while ladderset.iter().any(|&x| x) {
            ccs += 1;
            if ccs == b'Z' + 1 {
                ccs = b'a';
            }
            if ccs > b'z' {
                ccs = b'A';
            }
            let sheetset = self.findsheet(&mut ladderset);
            for i in 1..=nb {
                if sheetset[i] && self.bridges[i - 1].from == 0 {
                    if asci == 90 {
                        asci = 64;
                    }
                    asci += 1;
                    let name = if self.bridges[i - 1].btyp == BType::Parallel { asci + 32 } else { asci };
                    self.bridges[i - 1].laddername = name;
                    self.bridges[i - 1].sheetname = ccs;
                    self.bridges[i - 1].linkset = sheetset.clone();
                    let mut j = self.bridges[i - 1].towards;
                    while j != 0 {
                        self.bridges[j - 1].laddername = name;
                        self.bridges[j - 1].sheetname = ccs;
                        self.bridges[j - 1].linkset = sheetset.clone();
                        j = self.bridges[j - 1].towards;
                    }
                }
            }
        }
    }

    fn markstrands(&mut self) {
        let nb = self.bridges.len();
        for i in 1..=nb {
            if self.bridges[i - 1].from != 0 {
                continue;
            }
            // character sets of the beta1/beta2 columns on each strand
            let mut iset: [std::collections::BTreeSet<u8>; 2] = Default::default();
            let mut jset: [std::collections::BTreeSet<u8>; 2] = Default::default();
            let (mut ib0, mut ie0, mut jb0, mut je0) = (self.lchain, 0usize, self.lchain, 0usize);
            let mut j = i;
            loop {
                let w = self.bridges[j - 1].clone();
                for l in w.ib..=w.ie {
                    for b in 0..2 {
                        iset[b].insert(self.c[l].ss[BETA1 + b]);
                    }
                }
                for l in w.jb..=w.je {
                    for b in 0..2 {
                        jset[b].insert(self.c[l].ss[BETA1 + b]);
                    }
                }
                ib0 = ib0.min(w.ib);
                ie0 = ie0.max(w.ie);
                jb0 = jb0.min(w.jb);
                je0 = je0.max(w.je);
                j = w.towards;
                if j == 0 {
                    break;
                }
            }
            let blank: std::collections::BTreeSet<u8> = [b' '].into_iter().collect();
            // C tests column beta1 (index 0) for both choices
            let betai = if iset[0] == blank { 0 } else { 1 };
            let betaj = if jset[0] == blank { 0 } else { 1 };
            let mut j = i;
            loop {
                let w = self.bridges[j - 1].clone();
                for l in w.ib..=w.ie {
                    self.c[l].ss[BETA1 + betai] = w.laddername;
                    self.c[l].partner[betai] = if w.btyp == BType::Parallel {
                        (w.jb + l - w.ib) as i64
                    } else {
                        w.je as i64 - l as i64 + w.ib as i64
                    };
                }
                for l in w.jb..=w.je {
                    self.c[l].ss[BETA1 + betaj] = w.laddername;
                    self.c[l].partner[betaj] = if w.btyp == BType::Parallel {
                        (w.ib + l - w.jb) as i64
                    } else {
                        w.ie as i64 - l as i64 + w.jb as i64
                    };
                }
                j = w.towards;
                if j == 0 {
                    break;
                }
            }
            let cc = if ib0 == ie0 { b'B' } else { b'E' };
            for j in ib0..=ie0 {
                if self.c[j].ss[SYMBOL] != b'E' {
                    self.c[j].ss[SYMBOL] = cc;
                }
            }
            for j in jb0..=je0 {
                if self.c[j].ss[SYMBOL] != b'E' {
                    self.c[j].ss[SYMBOL] = cc;
                }
            }
        }
        for w in self.bridges.clone() {
            for l in w.ib..=w.ie {
                self.c[l].sheetlabel = w.sheetname;
            }
            for l in w.jb..=w.je {
                self.c[l].sheetlabel = w.sheetname;
            }
        }
    }

    fn flagbridge(&mut self) {
        self.bridges.clear();
        for i in 2..self.lchain {
            self.testbridge(i);
        }
        if self.bridges.is_empty() {
            return;
        }
        self.extendladder();
        self.sheet();
        self.markstrands();
    }

    fn flagsymbol(&mut self) {
        let nh = |c: u8| c == b'>' || c == b'X';
        let l = self.lchain as i64;
        let turn = |k: usize| TURN3 + k; // k: 0 -> turn3, 1 -> turn4, 2 -> turn5
        let mut i = 2i64;
        while i <= l - 4 {
            let iu = i as usize;
            if nh(self.c[iu - 1].ss[turn(1)]) && nh(self.c[iu].ss[turn(1)]) {
                for j in iu..=iu + 3 {
                    self.c[j].ss[SYMBOL] = b'H';
                }
            }
            i += 1;
        }
        let mut i = 2i64;
        while i <= l - 3 {
            let iu = i as usize;
            if nh(self.c[iu - 1].ss[turn(0)]) && nh(self.c[iu].ss[turn(0)]) {
                let empty = (iu..=iu + 2).all(|j| self.c[j].ss[SYMBOL] == b'G' || self.c[j].ss[SYMBOL] == b' ');
                if empty {
                    for j in iu..=iu + 2 {
                        self.c[j].ss[SYMBOL] = b'G';
                    }
                }
            }
            i += 1;
        }
        let mut i = 2i64;
        while i <= l - 5 {
            let iu = i as usize;
            if nh(self.c[iu - 1].ss[turn(2)]) && nh(self.c[iu].ss[turn(2)]) {
                let empty = (iu..=iu + 4).all(|j| self.c[j].ss[SYMBOL] == b'I' || self.c[j].ss[SYMBOL] == b' ');
                if empty {
                    for j in iu..=iu + 4 {
                        self.c[j].ss[SYMBOL] = b'I';
                    }
                }
            }
            i += 1;
        }
        for i in 2..self.lchain {
            if self.c[i].ss[SYMBOL] != b' ' {
                continue;
            }
            let mut cc = b' ';
            let mut j = 1;
            for t in 0..3 {
                j += 1;
                for k in 1..=j {
                    if i > k && nh(self.c[i - k].ss[turn(t)]) {
                        cc = b'T';
                    }
                }
            }
            if cc == b' ' {
                cc = self.c[i].ss[BEND];
            }
            self.c[i].ss[SYMBOL] = cc;
        }
    }

    fn flagturn(&mut self) {
        let mut k = 2usize;
        let mut cc = b'2';
        for t in 0..3 {
            k += 1;
            cc += 1;
            let col = TURN3 + t;
            let mut i = 1usize;
            while i + k <= self.lchain {
                if self.nochainbreak(i, i + k) && self.testbond(i + k, i) {
                    self.c[i + k].ss[col] = b'<';
                    for j in 1..k {
                        if self.c[i + j].ss[col] == b' ' {
                            self.c[i + j].ss[col] = cc;
                        }
                    }
                    self.c[i].ss[col] = if self.c[i].ss[col] == b'<' { b'X' } else { b'>' };
                }
                i += 1;
            }
        }
        for i in 1..=self.lchain {
            if self.c[i].kappa != 360.0 && self.c[i].kappa > 70.0 {
                self.c[i].ss[BEND] = b'S';
            }
        }
        self.flagsymbol();
    }
}

/// Run dsspcmbi on PDB text.
pub fn dssp_from_text(text: &str) -> Result<Dssp, DsspError> {
    let rd = input_coordinates(text)?;
    let mut st = State { c: rd.chain, lchain: rd.lchain, bridges: Vec::new() };
    st.flagchirality();
    st.flaghydrogenbonds();
    st.flagbridge();
    st.flagturn();
    let mut res = st.c;
    res.truncate(st.lchain + 1);
    Ok(Dssp { res })
}

/// Run dsspcmbi on a PDB file (.gz allowed).
pub fn dssp_from_pdb<P: AsRef<Path>>(path: P) -> Result<Dssp, DsspError> {
    let p = path.as_ref();
    let file = std::fs::File::open(p).map_err(|e| DsspError::Io(format!("{}: {}", p.display(), e)))?;
    let mut s = String::new();
    use std::io::Read;
    let r = if p.extension().map_or(false, |x| x == "gz") {
        flate2::read::GzDecoder::new(file).read_to_string(&mut s)
    } else {
        let mut f = file;
        f.read_to_string(&mut s)
    };
    r.map_err(|e| DsspError::Io(format!("{}: {}", p.display(), e)))?;
    dssp_from_text(&s)
}

impl Dssp {
    /// What util.f getdssp hands puu for one chain: summary structure per residue
    /// (breaks skipped) and bridge partners renumbered to that sequence, -99 for a
    /// partner in another chain.
    pub fn to_puu_input(&self, chainid: u8) -> super::puu::DsspInput {
        let mut fssp = vec![0i32; self.res.len()];
        let mut struc = Vec::new();
        let mut raw = Vec::new();
        for i in 1..self.res.len() {
            let r = &self.res[i];
            if r.is_break() {
                continue;
            }
            if chainid != b'*' && r.chain_id() != chainid {
                continue;
            }
            struc.push(r.ss[SYMBOL]);
            raw.push(r.partner);
            fssp[i] = struc.len() as i32;
        }
        let bp = raw
            .iter()
            .map(|p| {
                let mut o = [0i32; 2];
                for k in 0..2 {
                    o[k] = if p[k] > 0 {
                        match fssp.get(p[k] as usize).copied().unwrap_or(0) {
                            x if x > 0 => x,
                            _ => -99,
                        }
                    } else {
                        p[k] as i32
                    };
                }
                o
            })
            .collect();
        super::puu::DsspInput { struc, bp }
    }
}
