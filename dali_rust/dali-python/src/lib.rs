use std::sync::Arc;

use numpy::PyArray2;
use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;

use dali_core::io::{dat, import};
use dali_core::numerics::compute_transform;
use dali_core::pipeline;
use dali_core::store::ProteinStore;
use dali_core::{AlignmentBlock, DccpEntry, Protein};

// ── PyAlignmentBlock ──────────────────────────────────────────────

#[pyclass(name = "AlignmentBlock")]
#[derive(Clone)]
struct PyAlignmentBlock {
    #[pyo3(get)]
    l1: u32,
    #[pyo3(get)]
    r1: u32,
    #[pyo3(get)]
    l2: u32,
    #[pyo3(get)]
    r2: u32,
}

#[pymethods]
impl PyAlignmentBlock {
    fn __repr__(&self) -> String {
        format!("AlignmentBlock(l1={}, r1={}, l2={}, r2={})",
                self.l1, self.r1, self.l2, self.r2)
    }
}

impl From<&AlignmentBlock> for PyAlignmentBlock {
    fn from(b: &AlignmentBlock) -> Self {
        PyAlignmentBlock { l1: b.l1, r1: b.r1, l2: b.l2, r2: b.r2 }
    }
}

// ── PyDccpEntry ───────────────────────────────────────────────────

#[pyclass(name = "DccpEntry")]
struct PyDccpEntry {
    inner: DccpEntry,
}

#[pymethods]
impl PyDccpEntry {
    #[getter]
    fn cd1(&self) -> &str { &self.inner.cd1 }

    #[getter]
    fn cd2(&self) -> &str { &self.inner.cd2 }

    #[getter]
    fn score(&self) -> f64 { self.inner.score }

    #[getter]
    fn zscore(&self) -> f64 { self.inner.zscore }

    #[getter]
    fn rmsd(&self) -> f64 { self.inner.rmsd }

    #[getter]
    fn nblock(&self) -> usize { self.inner.blocks.len() }

    #[getter]
    fn blocks(&self) -> Vec<PyAlignmentBlock> {
        self.inner.blocks.iter().map(PyAlignmentBlock::from).collect()
    }

    fn __repr__(&self) -> String {
        format!("DccpEntry(cd1='{}', cd2='{}', score={:.1}, zscore={:.1}, rmsd={:.1}, nblock={})",
                self.inner.cd1, self.inner.cd2, self.inner.score,
                self.inner.zscore, self.inner.rmsd, self.inner.blocks.len())
    }
}

// ── PyProtein ─────────────────────────────────────────────────────

#[pyclass(name = "Protein")]
struct PyProtein {
    inner: Arc<Protein>,
}

#[pymethods]
impl PyProtein {
    #[getter]
    fn code(&self) -> &str { &self.inner.code }

    #[getter]
    fn nres(&self) -> usize { self.inner.nres }

    #[getter]
    fn nseg(&self) -> usize { self.inner.nseg }

    #[getter]
    fn na(&self) -> usize { self.inner.na }

    #[getter]
    fn nb(&self) -> usize { self.inner.nb }

    #[getter]
    fn sequence(&self) -> &str { &self.inner.sequence }

    /// Return CA coordinates as a numpy array of shape (3, nres).
    fn ca_coords<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        PyArray2::from_array(py, &self.inner.ca)
    }

    /// Return PDB residue serial numbers (1-based).
    /// For .dat-loaded proteins, this is sequential 1..=nres.
    /// For PDB-imported proteins, this reflects actual PDB numbering.
    #[getter]
    fn resid_map(&self) -> Vec<i32> {
        self.inner.resid_map.clone()
    }

    /// Write this protein to a .dat file.
    fn write_dat(&self, path: &str) -> PyResult<()> {
        dat::write_dat(&self.inner, path)
            .map_err(|e| PyValueError::new_err(format!("{}", e)))
    }

    fn __repr__(&self) -> String {
        format!("Protein(code='{}', nres={}, nseg={})",
                self.inner.code, self.inner.nres, self.inner.nseg)
    }
}

// ── PyProteinStore ────────────────────────────────────────────────

#[pyclass(name = "ProteinStore")]
struct PyProteinStore {
    inner: ProteinStore,
}

#[pymethods]
impl PyProteinStore {
    #[new]
    fn new(dat_dir: &str) -> Self {
        PyProteinStore { inner: ProteinStore::new(dat_dir) }
    }

    fn get_protein(&self, code: &str) -> PyResult<PyProtein> {
        self.inner.get_protein(code)
            .map(|p| PyProtein { inner: p })
            .map_err(|e| PyValueError::new_err(format!("{:?}", e)))
    }

    fn __len__(&self) -> usize { self.inner.len() }

    fn __contains__(&self, code: &str) -> bool {
        // Check filesystem, not just cache
        let path = std::path::Path::new(self.inner.dat_dir()).join(format!("{}.dat", code));
        path.exists()
    }

    fn __repr__(&self) -> String {
        format!("ProteinStore(dat_dir='{}', loaded={})",
                self.inner.dat_dir(), self.inner.len())
    }
}

// ── Module-level functions ────────────────────────────────────────

/// Compare two proteins through the full DALI pipeline.
#[pyfunction]
fn compare_pair(cd1: &str, cd2: &str, store: &PyProteinStore) -> Vec<PyDccpEntry> {
    pipeline::compare_pair(cd1, cd2, &store.inner)
        .into_iter()
        .map(|e| PyDccpEntry { inner: e })
        .collect()
}

/// Run the WOLF→DP→DALICON path for a query against multiple targets.
#[pyfunction]
fn run_wolf_path(query: &str, targets: Vec<String>, store: &PyProteinStore) -> Vec<PyDccpEntry> {
    let refs: Vec<&str> = targets.iter().map(|s| s.as_str()).collect();
    pipeline::run_wolf_path(query, &refs, &store.inner)
        .into_iter()
        .map(|e| PyDccpEntry { inner: e })
        .collect()
}

/// Run the PARSI→FILTER95→DALICON→DP path for a query against multiple targets.
#[pyfunction]
fn run_parsi_path(query: &str, targets: Vec<String>, store: &PyProteinStore) -> Vec<PyDccpEntry> {
    let refs: Vec<&str> = targets.iter().map(|s| s.as_str()).collect();
    pipeline::run_parsi_path(query, &refs, &store.inner)
        .into_iter()
        .map(|e| PyDccpEntry { inner: e })
        .collect()
}

/// Read a .dat file and return a Protein object.
#[pyfunction]
fn read_dat(path: &str) -> PyResult<PyProtein> {
    dat::read_dat(path)
        .map(|p| PyProtein { inner: Arc::new(p) })
        .map_err(|e| PyValueError::new_err(format!("{:?}", e)))
}

/// Import a PDB/CIF file and return a Protein object.
#[pyfunction]
fn import_pdb(path: &str, chain: &str, pdb_code: &str) -> PyResult<PyProtein> {
    import::import_pdb(path, chain, pdb_code)
        .map(|p| PyProtein { inner: Arc::new(p) })
        .map_err(|e| PyValueError::new_err(format!("{:?}", e)))
}

// ── PyAlignResult ────────────────────────────────────────────────

/// Result of a PDB-to-PDB structural alignment.
#[pyclass(name = "AlignResult")]
struct PyAlignResult {
    #[pyo3(get)]
    zscore: f64,
    #[pyo3(get)]
    score: f64,
    #[pyo3(get)]
    rmsd: f64,
    #[pyo3(get)]
    n_aligned: usize,
    #[pyo3(get)]
    alignments: Vec<(i32, i32)>,  // (query_pdb_resid, template_pdb_resid)
    #[pyo3(get)]
    rotation: Vec<Vec<f64>>,       // 3x3 rotation matrix
    #[pyo3(get)]
    translation: Vec<f64>,         // 3-element translation vector
    #[pyo3(get)]
    blocks: Vec<PyAlignmentBlock>,
}

#[pymethods]
impl PyAlignResult {
    fn __repr__(&self) -> String {
        format!("AlignResult(zscore={:.1}, score={:.1}, rmsd={:.1}, n_aligned={})",
                self.zscore, self.score, self.rmsd, self.n_aligned)
    }
}

/// Align two PDB/CIF files and return the best structural alignment.
///
/// Imports both structures, runs the full DALI pipeline, and returns
/// the best hit (highest Z-score) with rotation, translation, and
/// residue-level alignment pairs in PDB numbering.
#[pyfunction]
#[pyo3(signature = (query_path, template_path, query_chain="A", template_chain="A", query_code="query", template_code="templ"))]
fn align_pdb(
    query_path: &str,
    template_path: &str,
    query_chain: &str,
    template_chain: &str,
    query_code: &str,
    template_code: &str,
) -> PyResult<Option<PyAlignResult>> {
    // Import both proteins
    let query = import::import_pdb(query_path, query_chain, query_code)
        .map_err(|e| PyValueError::new_err(format!("Query import error: {}", e)))?;
    let template = import::import_pdb(template_path, template_chain, template_code)
        .map_err(|e| PyValueError::new_err(format!("Template import error: {}", e)))?;

    // Write to temp .dat files
    let tmp_dir = std::env::temp_dir().join(format!("dali_align_{}", std::process::id()));
    std::fs::create_dir_all(&tmp_dir)
        .map_err(|e| PyValueError::new_err(format!("Cannot create temp dir: {}", e)))?;

    let q_dat = tmp_dir.join(format!("{}.dat", query.code));
    let t_dat = tmp_dir.join(format!("{}.dat", template.code));

    let cleanup = || { let _ = std::fs::remove_dir_all(&tmp_dir); };

    dat::write_dat(&query, &q_dat)
        .map_err(|e| { cleanup(); PyValueError::new_err(format!("Write query .dat: {}", e)) })?;
    dat::write_dat(&template, &t_dat)
        .map_err(|e| { cleanup(); PyValueError::new_err(format!("Write template .dat: {}", e)) })?;

    // Run pipeline
    let store = ProteinStore::new(&tmp_dir);
    let results = pipeline::compare_pair(&query.code, &template.code, &store);

    // Filter to query→template direction only (blocks l1/r1 refer to query)
    let q_results: Vec<&DccpEntry> = results.iter()
        .filter(|r| r.cd1 == query.code && r.cd2 == template.code)
        .collect();

    if q_results.is_empty() {
        cleanup();
        return Ok(None);
    }

    // Pick best result by Z-score
    let best = q_results.iter()
        .max_by(|a, b| a.zscore.partial_cmp(&b.zscore).unwrap())
        .unwrap();

    // Compute rotation + translation from alignment blocks
    // blocks l1/r1 = query indices, l2/r2 = template indices
    let (rotation, translation) = match compute_transform(&query.ca, &template.ca, &best.blocks) {
        Some((u, t)) => {
            let rot: Vec<Vec<f64>> = (0..3).map(|i| {
                (0..3).map(|j| u[[i, j]]).collect()
            }).collect();
            let trans: Vec<f64> = t.to_vec();
            (rot, trans)
        }
        None => {
            cleanup();
            return Ok(None);
        }
    };

    // Expand blocks to residue-level pairs in PDB numbering
    let mut alignments: Vec<(i32, i32)> = Vec::new();
    let mut n_aligned = 0usize;
    for b in &best.blocks {
        let len = (b.r1 - b.l1 + 1) as usize;
        for k in 0..len {
            let q_idx = (b.l1 as usize - 1) + k;  // 0-based internal index
            let t_idx = (b.l2 as usize - 1) + k;
            if q_idx < query.resid_map.len() && t_idx < template.resid_map.len() {
                alignments.push((query.resid_map[q_idx], template.resid_map[t_idx]));
                n_aligned += 1;
            }
        }
    }

    let blocks: Vec<PyAlignmentBlock> = best.blocks.iter().map(PyAlignmentBlock::from).collect();

    cleanup();

    Ok(Some(PyAlignResult {
        zscore: best.zscore,
        score: best.score,
        rmsd: best.rmsd,
        n_aligned,
        alignments,
        rotation,
        translation,
        blocks,
    }))
}

/// Write a Protein to a .dat file (standalone function).
#[pyfunction]
fn write_dat(protein: &PyProtein, path: &str) -> PyResult<()> {
    dat::write_dat(&protein.inner, path)
        .map_err(|e| PyValueError::new_err(format!("{}", e)))
}

// ── Module definition ─────────────────────────────────────────────

#[pymodule]
fn dali(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyAlignmentBlock>()?;
    m.add_class::<PyDccpEntry>()?;
    m.add_class::<PyProtein>()?;
    m.add_class::<PyProteinStore>()?;
    m.add_class::<PyAlignResult>()?;
    m.add_function(wrap_pyfunction!(compare_pair, m)?)?;
    m.add_function(wrap_pyfunction!(run_wolf_path, m)?)?;
    m.add_function(wrap_pyfunction!(run_parsi_path, m)?)?;
    m.add_function(wrap_pyfunction!(read_dat, m)?)?;
    m.add_function(wrap_pyfunction!(import_pdb, m)?)?;
    m.add_function(wrap_pyfunction!(align_pdb, m)?)?;
    m.add_function(wrap_pyfunction!(write_dat, m)?)?;
    Ok(())
}
