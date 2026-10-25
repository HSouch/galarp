
//! TOML Configuration files for GalaRP
//! 
//! Reading, writing, and validation


use std::fmt;
use std::path::{Path, PathBuf};

use bigwig::cluster::{hse, BetaModel, ClusterGas, RadialTable};
use bigwig::orbit::{self, Body, Motion, OrbitTable};
use bigwig::potentials::{Burkert, Composite, Hernquist, MiyamotoNagai, Nfw, Param, PointMass, Potential};
use bigwig::units::{self, UNIT_DENSITY_CGS, UNIT_TIME_S};
use serde::Deserialize;


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


#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Config {
    pub run: RunConfig,
   
}


#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RunConfig {
    /// Where snapshots, the log and a copy of the config go. Relative paths are relative to
    /// the directory the program is run from.
    pub output_dir: PathBuf,
    pub t_end_myr: f64,
    /// Log interval; a whole multiple of `restricted.dt_myr`.
    #[serde(default = "d_dt_max")]
    pub dt_max_myr: f64,
    /// Must be a whole multiple of `dt_max_myr`.
    #[serde(default = "d_snapshot")]
    pub snapshot_interval_myr: f64,
    /// Gzip + shuffle on snapshot datasets.
    #[serde(default)]
    pub compress: bool,
    /// Mean molecular weight written to snapshot headers.
    #[serde(default = "d_mu")]
    pub mean_molecular_weight: f64,
}

fn d_dt_max() -> f64 {
    2.0
}
fn d_snapshot() -> f64 {
    10.0
}
fn d_mu() -> f64 {
    0.6
}