use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::{Arc, RwLock};

use crate::io::dat;
use crate::types::{DistMatScale10, DistMatScale100, Protein};
use crate::numerics::scoring::dpgetdist;

/// Thread-safe protein data cache.
///
/// Holds immutable data behind Arc. Modules that need to mutate coordinates
/// (e.g. DALICON's ca1, FILTER95's xca) must clone from the store.
pub struct ProteinStore {
    dat_dir: PathBuf,
    proteins: RwLock<HashMap<String, Arc<Protein>>>,
    dist_scale10: RwLock<HashMap<String, Arc<DistMatScale10>>>,
    dist_scale100: RwLock<HashMap<String, Arc<DistMatScale100>>>,
}

impl ProteinStore {
    /// Create a new store that reads .dat files from the given directory.
    pub fn new<P: AsRef<Path>>(dat_dir: P) -> Self {
        ProteinStore {
            dat_dir: dat_dir.as_ref().to_path_buf(),
            proteins: RwLock::new(HashMap::new()),
            dist_scale10: RwLock::new(HashMap::new()),
            dist_scale100: RwLock::new(HashMap::new()),
        }
    }

    /// Get or load a protein by code (e.g. "101mA").
    pub fn get_protein(&self, code: &str) -> Result<Arc<Protein>, dat::DatError> {
        // Check cache first (read lock)
        {
            let cache = self.proteins.read().unwrap();
            if let Some(prot) = cache.get(code) {
                return Ok(Arc::clone(prot));
            }
        }

        // Load from file (write lock)
        let filepath = self.dat_dir.join(format!("{}.dat", code));
        let protein = dat::read_dat(&filepath)?;
        let arc = Arc::new(protein);

        let mut cache = self.proteins.write().unwrap();
        cache.insert(code.to_string(), Arc::clone(&arc));
        Ok(arc)
    }

    /// Get or compute distance matrix at scale 10 (for DP module).
    pub fn get_dist_scale10(&self, code: &str) -> Result<Arc<DistMatScale10>, dat::DatError> {
        // Check cache
        {
            let cache = self.dist_scale10.read().unwrap();
            if let Some(d) = cache.get(code) {
                return Ok(Arc::clone(d));
            }
        }

        // Compute from protein
        let prot = self.get_protein(code)?;
        let data = dpgetdist(&prot.ca, prot.nres);
        let dist = DistMatScale10 {
            data,
            nres: prot.nres,
        };
        let arc = Arc::new(dist);

        let mut cache = self.dist_scale10.write().unwrap();
        cache.insert(code.to_string(), Arc::clone(&arc));
        Ok(arc)
    }

    /// Get or compute distance matrix at scale 100 (for FILTER95 module).
    pub fn get_dist_scale100(&self, code: &str) -> Result<Arc<DistMatScale100>, dat::DatError> {
        // Check cache
        {
            let cache = self.dist_scale100.read().unwrap();
            if let Some(d) = cache.get(code) {
                return Ok(Arc::clone(d));
            }
        }

        // Compute from protein (scale 100 instead of 10)
        let prot = self.get_protein(code)?;
        let nres = prot.nres;
        let mut data = ndarray::Array2::<i16>::zeros((nres, nres));
        let hundred: f32 = 100.0;

        for i in 0..nres {
            let xi = prot.ca[[0, i]] as f32;
            let yi = prot.ca[[1, i]] as f32;
            let zi = prot.ca[[2, i]] as f32;
            for j in 0..i {
                let dx = xi - prot.ca[[0, j]] as f32;
                let dy = yi - prot.ca[[1, j]] as f32;
                let dz = zi - prot.ca[[2, j]] as f32;
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                let x = crate::numerics::nint_f32(hundred * dist).min(32767);
                data[[i, j]] = x as i16;
                data[[j, i]] = x as i16;
            }
        }

        let dist = DistMatScale100 { data, nres };
        let arc = Arc::new(dist);

        let mut cache = self.dist_scale100.write().unwrap();
        cache.insert(code.to_string(), Arc::clone(&arc));
        Ok(arc)
    }

    /// Check if a protein is already cached.
    pub fn has_protein(&self, code: &str) -> bool {
        self.proteins.read().unwrap().contains_key(code)
    }

    /// Number of cached proteins.
    pub fn len(&self) -> usize {
        self.proteins.read().unwrap().len()
    }

    /// Is the store empty?
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Get the dat directory path as a string.
    pub fn dat_dir(&self) -> &str {
        self.dat_dir.to_str().unwrap_or("")
    }
}
