//! `[cluster]`: the cluster's potential and its gas (the ICM).

use std::path::{Path, PathBuf};

use galarp::cluster::{hse, BetaModel, ClusterGas, RadialTable};
use galarp::potentials::Composite;
use galarp::units::{self, UNIT_DENSITY_CGS};
use serde::Deserialize;

use super::{composite, d_mu, err, ConfigError, PotentialTerm};

// ---------------------------------------------------------------------------
// Schema
// ---------------------------------------------------------------------------

/// `[cluster]`
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ClusterConfig {
    pub gas: ClusterGasConfig,
    /// Cluster potential terms (empty: no cluster gravity, as for a uniform wind).
    #[serde(default)]
    pub potential: Vec<PotentialTerm>,
}

/// `[cluster.gas]`, picked by `model`.
#[derive(Clone, Debug, Deserialize)]
#[serde(tag = "model", rename_all = "snake_case", deny_unknown_fields)]
pub enum ClusterGasConfig {
    /// No ambient gas at all: an isolated galaxy in vacuum (no ram pressure possible).
    None,
    Uniform {
        density_g_cm3: f64,
        temperature_k: f64,
        #[serde(default = "d_mu")]
        mean_molecular_weight: f64,
    },
    /// Beta model in hydrostatic equilibrium with the cluster potential.
    Beta {
        rho0_g_cm3: f64,
        r_c_kpc: f64,
        beta: f64,
        #[serde(default = "d_table_rmin")]
        r_min_kpc: f64,
        #[serde(default = "d_table_rmax")]
        r_max_kpc: f64,
        #[serde(default = "d_table_n")]
        n_table: usize,
    },
    /// Radial profile from a CSV with columns `r_kpc,density_g_cm3,temperature_k`. With
    /// `hydrostatic = true`, the temperature column is ignored and recomputed from
    /// hydrostatic equilibrium in the cluster potential.
    Table {
        file: PathBuf,
        #[serde(default = "d_mu")]
        mean_molecular_weight: f64,
        #[serde(default)]
        hydrostatic: bool,
    },
}

// ---------------------------------------------------------------------------
// Defaults
// ---------------------------------------------------------------------------

fn d_table_rmin() -> f64 {
    1.0
}

fn d_table_rmax() -> f64 {
    1e5
}

fn d_table_n() -> usize {
    4000
}

/// Radial points in the hydrostatic table built from a `table` profile.
const HSE_TABLE_N: usize = 4000;

// ---------------------------------------------------------------------------
// Validation (this section only; whether a potential is required is checked in
// `Config::validate`, since integrated orbits need one too)
// ---------------------------------------------------------------------------

impl ClusterConfig {
    pub(super) fn validate(&self) -> Result<(), ConfigError> {
        for t in &self.potential {
            t.validate()?;
        }
        Ok(())
    }
}

impl ClusterGasConfig {
    /// Gas models that are built in hydrostatic equilibrium with the cluster potential.
    pub fn needs_potential(&self) -> bool {
        matches!(self, ClusterGasConfig::Beta { .. } | ClusterGasConfig::Table { hydrostatic: true, .. })
    }
}

// ---------------------------------------------------------------------------
// Conversion to code units / galarp types
// ---------------------------------------------------------------------------

impl ClusterConfig {
    /// The cluster potential, or `None` if no terms are given.
    pub fn build_potential(&self) -> Option<Composite> {
        (!self.potential.is_empty()).then(|| composite(&self.potential))
    }

    /// The ICM in code units, for adiabatic index `gamma`.
    pub fn build_gas(&self, gamma: f64) -> Result<ClusterGas, ConfigError> {
        match &self.gas {
            ClusterGasConfig::None => Ok(ClusterGas::Uniform { rho: 0.0, u: 0.0 }),
            ClusterGasConfig::Uniform { density_g_cm3, temperature_k, mean_molecular_weight } => Ok(ClusterGas::Uniform {
                rho: density_g_cm3 / UNIT_DENSITY_CGS,
                u: units::u_from_temperature(*temperature_k, *mean_molecular_weight, gamma),
            }),
            ClusterGasConfig::Beta { rho0_g_cm3, r_c_kpc, beta, r_min_kpc, r_max_kpc, n_table } => {
                let model = BetaModel { rho_0: rho0_g_cm3 / UNIT_DENSITY_CGS, r_c: *r_c_kpc, beta: *beta };
                let pot = self.build_potential().expect("validated");
                Ok(ClusterGas::beta_model_hse(&model, &pot, gamma, *r_min_kpc, *r_max_kpc, *n_table))
            }
            ClusterGasConfig::Table { file, mean_molecular_weight, hydrostatic } => {
                let (r, rho, t) = read_profile(file)?;
                let rho: Vec<f64> = rho.iter().map(|v| v / UNIT_DENSITY_CGS).collect();
                if *hydrostatic {
                    let pot = self.build_potential().expect("validated");
                    let density = RadialTable::new(&r, &rho, &vec![1.0; r.len()]);
                    let (r_min, r_max) = (r[0], r[r.len() - 1]);
                    let table = hse::hydrostatic_table(|x| density.density(x), &pot, gamma, r_min, r_max, HSE_TABLE_N);
                    Ok(ClusterGas::Radial(table))
                } else {
                    Ok(ClusterGas::Radial(RadialTable::from_temperature(&r, &rho, &t, *mean_molecular_weight, gamma)))
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Profile CSV
// ---------------------------------------------------------------------------

const PROFILE_HEADER: &str = "r_kpc,density_g_cm3,temperature_k";

/// Reads a `r_kpc,density_g_cm3,temperature_k` CSV (header line required; `#` comments).
fn read_profile(path: &Path) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>), ConfigError> {
    let text = std::fs::read_to_string(path).map_err(|e| ConfigError(format!("{}: {e}", path.display())))?;
    let (mut r, mut rho, mut t) = (Vec::new(), Vec::new(), Vec::new());
    let mut header = false;
    for (n, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if !header {
            if line.replace(' ', "") != PROFILE_HEADER {
                return err(format!("{}: expected header {PROFILE_HEADER}", path.display()));
            }
            header = true;
            continue;
        }
        let vals: Result<Vec<f64>, _> = line.split(',').map(|s| s.trim().parse::<f64>()).collect();
        match vals {
            Ok(v) if v.len() == 3 => {
                r.push(v[0]);
                rho.push(v[1]);
                t.push(v[2]);
            }
            _ => return err(format!("{} line {}: expected three numbers", path.display(), n + 1)),
        }
    }
    if r.len() < 2 {
        return err(format!("{}: need at least two rows", path.display()));
    }
    Ok((r, rho, t))
}
