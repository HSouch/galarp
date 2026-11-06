//! TOML configuration files for GalaRP: 
//! Reading and validation.
//!
//! Each section of the config file has its own file.
//!
//! | file              | sections                                                 |
//! |-------------------|----------------------------------------------------------|
//! | `run.rs`          | `[run]`, `[hydro]`, `[box]`, `[galarp]`, `[diagnostics]` |
//! | `cluster.rs`      | `[cluster]`, `[cluster.gas]`                             |
//! | `galaxy.rs`       | `[[galaxy]]`, `[galaxy.disk]`                            |
//! | `orbit.rs`        | `[galaxy.orbit]`, `[orbits]`                             |
//! | `potential.rs`    | `[[galaxy.potential]]`, `[[cluster.potential]]`          |
//! | `param.rs`        | `ParamSpec`                                              |
//! | `ram_pressure.rs` | `[ram_pressure]`, `[ram_pressure.shadow]`                |
//!
//! 
//! To add a parameter:
//! 1. add the field (named with unit suffix _UNIT)
//! 2. Choose whether required or optional. A default is generally preferable.
//!     #[serde(default = "d_DEFAULT_FN")]
//! 3. If the file can change over time, make it a ParamSpec
//! 4. Add validation to the struct's validate.
//! 
//! 



use std::fmt;
use std::path::{Path, PathBuf};

use galarp::cluster::ClusterGas;
use galarp::potentials::Composite;
use serde::Deserialize;



mod galaxy;
pub use galaxy::{DiskConfig, GalaxyConfig};

mod orbit;
pub use orbit::{OrbitConfig, OrbitsConfig};

mod param;
pub use param::ParamSpec;

mod potential;
pub use potential::PotentialTerm;
use potential::composite;

mod cluster;
pub use cluster::{ClusterConfig, ClusterGasConfig};

mod ram_pressure;
pub use ram_pressure::{RamPressureConfig, ShadowConfig};

mod run;
pub use run::{BoxSection, DiagnosticsConfig, GalaRPConfig, HydroConfig, RunConfig};
use run::d_mu;

#[cfg(test)]
mod tests;


#[derive(Debug)]
pub struct ConfigError(pub String);

impl fmt::Display for ConfigError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for ConfigError {}

fn err<T>(msg: impl Into<String>) -> Result<T, ConfigError> {
    Err(ConfigError(msg.into()))
}

impl ConfigError {
    /// Prefixes the message with where the error is, e.g. "galaxy 1 (companion)".
    fn at(self, prefix: impl fmt::Display) -> Self {
        ConfigError(format!("{prefix}: {}", self.0))
    }
}

/// `Ok` if `v` is positive and finite.
fn positive(name: &str, v: f64) -> Result<(), ConfigError> {
    if v > 0.0 && v.is_finite() { Ok(()) } else { err(format!("{name} must be positive, got {v}")) }
}



#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Config {
    pub run: RunConfig,

    #[serde(default)]
    pub hydro: HydroConfig,

    #[serde(rename = "box")]
    pub box_: BoxSection,
    pub cluster: ClusterConfig,

    pub galaxy: Vec<GalaxyConfig>,

    #[serde(default)]
    pub diagnostics: DiagnosticsConfig,

    #[serde(default)]
    pub orbits: OrbitsConfig,

    /// Analytic ram-pressure wind config on the clouds.
    #[serde(default)]
    pub ram_pressure: RamPressureConfig,

    #[serde(default)]
    pub galarp: GalaRPConfig,
    
    /// A copy of the run's state appended to each output's  config. 
    /// Ignored on input, so such a copy can still be read.
    #[serde(default)]
    pub output: Option<toml::Table>,
}


impl Config {
    pub fn from_str_in(text: &str, base_dir: &Path) -> Result<Self, ConfigError> {
        let mut cfg: Config = toml::from_str(text).map_err(|e| ConfigError(format!("config: {e}")))?;
        let resolve = |p: &mut PathBuf| {
            if p.is_relative() {
                *p = base_dir.join(&*p);
            }
        };
        if let ClusterGasConfig::Table { file, .. } = &mut cfg.cluster.gas {
            resolve(file);
        }
        for g in &mut cfg.galaxy {
            if let OrbitConfig::Table { file } = &mut g.orbit {
                resolve(file);
            }
        }
        cfg.validate()?;
        Ok(cfg)
    }

    pub fn from_file(path: impl AsRef<Path>) -> Result<Self, ConfigError> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path).map_err(|e| ConfigError(format!("{}: {e}", path.display())))?;
        let base = path.parent().unwrap_or(Path::new("."));
        Self::from_str_in(&text, base).map_err(|e| ConfigError(format!("{}: {}", path.display(), e.0)))
    }

    fn validate(&self) -> Result<(), ConfigError> {
        self.run.validate()?;
        self.box_.validate()?;
        self.hydro.validate()?;

        if self.galaxy.is_empty() {
            return err("need at least one [[galaxy]] in");
        }

        for (k, g) in self.galaxy.iter().enumerate() {
            g.validate().map_err(|e| e.at(format!("galaxy {k} ({})", g.label(k))))?;
        }

        self.cluster.validate()?;
        self.ram_pressure.validate()?;
        self.galarp.validate()?;

        let pulled_by_galaxies = self.orbits.mutual_gravity && self.galaxy.len() > 1;
        let needs_cluster_potential = self.cluster.gas.needs_potential()
            || (!pulled_by_galaxies && self.galaxy.iter().any(|g| g.orbit.is_integrated()));
        if needs_cluster_potential && self.cluster.potential.is_empty() {
            return err("this cluster gas model or orbit mode needs [[cluster.potential]] terms");
        }

        if self.ram_pressure.enabled && matches!(self.cluster.gas, ClusterGasConfig::None) {
            return err("ram_pressure needs ICM gas: [cluster.gas]")
        }

        Ok(())
    }
}


impl Config {

    /// Label of galaxy `k`: its `name`, or "galaxy k".
    pub fn galaxy_name(&self, k: usize) -> String {
        self.galaxy[k].label(k)
    }

    /// Initial-condition options for galaxy `k`'s disk (self-gravity, noise).
    pub fn disk_options(&self, k: usize) -> galarp::ic::DiskOptions {
        self.galaxy[k].disk_options()
    }

    /// Orientation of galaxy `k` (its disk and potential).
    pub fn galaxy_orientation(&self, k: usize) -> galarp::orientation::Orientation {
        self.galaxy[k].orientation()
    }

    /// Disk parameters (code units) and cloud count of galaxy `k`, if it has a disk.
    pub fn disk(&self, k: usize) -> Option<(galarp::ic::GasDiskParams, usize)> {
        self.galaxy[k].disk.as_ref().map(|d| (d.params(), d.n_particles))
    }

    /// Cluster potential, or `None` if no terms are given.
    pub fn cluster_potential(&self) -> Option<Composite> {
        self.cluster.build_potential()
    }

    /// Potential of galaxy `k` (0 = parent).
    pub fn galaxy_potential(&self, k: usize) -> Composite {
        self.galaxy[k].build_potential()
    }

    pub fn cluster_gas(&self) -> Result<ClusterGas, ConfigError> {
        self.cluster.build_gas(self.hydro.gamma)
    }
}