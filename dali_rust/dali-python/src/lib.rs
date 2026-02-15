use std::sync::Arc;

use numpy::PyArray2;
use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;

use dali_core::io::{dat, import};
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

// ── Module definition ─────────────────────────────────────────────

#[pymodule]
fn dali(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyAlignmentBlock>()?;
    m.add_class::<PyDccpEntry>()?;
    m.add_class::<PyProtein>()?;
    m.add_class::<PyProteinStore>()?;
    m.add_function(wrap_pyfunction!(compare_pair, m)?)?;
    m.add_function(wrap_pyfunction!(run_wolf_path, m)?)?;
    m.add_function(wrap_pyfunction!(run_parsi_path, m)?)?;
    m.add_function(wrap_pyfunction!(read_dat, m)?)?;
    m.add_function(wrap_pyfunction!(import_pdb, m)?)?;
    Ok(())
}
