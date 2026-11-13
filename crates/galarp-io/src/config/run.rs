//! `[run]`, `[hydro]`, `[box]`, `[galarp]` and `[diagnostics]`: run length, time steps
//! and other settings that don't belong to the cluster or a galaxy.

use std::path::PathBuf;

use galarp::units::{self, MYR};
use serde::Deserialize;

use super::{err, positive, Config, ConfigError};

// ---------------------------------------------------------------------------
// Schema
// ---------------------------------------------------------------------------

/// `[run]`
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RunConfig {
    /// Output directory. Relative paths are from the directory the program is run from.
    pub output_dir: PathBuf,
    pub t_end_myr: f64,
    /// Log interval.
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

/// `[hydro]`
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct HydroConfig {
    pub gamma: f64,
}

/// `[box]`
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BoxSection {
    pub side_kpc: f64,
}

/// `[galarp]`: GalaRP-specific integration settings.
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct GalaRPConfig {
    /// RK4 step; must divide `run.dt_max_myr` (the log interval).
    pub dt_myr: f64,
    /// Vertical velocity dispersion of the clouds. Default: each disk's `sound_speed_km_s`,
    /// which keeps the disk's vertical profile.
    pub vertical_dispersion_km_s: Option<f64>,
    /// Seed of the vertical velocities.
    pub seed: u64,
}

/// `[diagnostics]`: the region counted as "still in the disk" for the log's
/// `in_disk_fraction` column, a cylinder around the galaxy center aligned with the disk.
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct DiagnosticsConfig {
    pub disk_region_radius_kpc: f64,
    pub disk_region_half_height_kpc: f64,
}

// ---------------------------------------------------------------------------
// Defaults
// ---------------------------------------------------------------------------

fn d_dt_max() -> f64 {
    2.0
}

fn d_snapshot() -> f64 {
    10.0
}

/// Mean molecular weight; also the default for the cluster gas models.
pub(super) fn d_mu() -> f64 {
    0.6
}

impl Default for HydroConfig {
    fn default() -> Self {
        Self { gamma: units::GAMMA }
    }
}

impl Default for GalaRPConfig {
    fn default() -> Self {
        Self { dt_myr: 0.5, vertical_dispersion_km_s: None, seed: 1 }
    }
}

impl Default for DiagnosticsConfig {
    fn default() -> Self {
        Self { disk_region_radius_kpc: 15.0, disk_region_half_height_kpc: 2.0 }
    }
}

// ---------------------------------------------------------------------------
// Validation
// ---------------------------------------------------------------------------

/// `a / b` as a whole number of at least 1, or `None` if it isn't one.
fn whole_ratio(a: f64, b: f64) -> Option<u64> {
    let ratio = a / b;
    ((ratio - ratio.round()).abs() <= 1e-9 && ratio.round() >= 1.0).then(|| ratio.round() as u64)
}

impl RunConfig {
    pub(super) fn validate(&self) -> Result<(), ConfigError> {
        positive("run.t_end_myr", self.t_end_myr)?;
        positive("run.dt_max_myr", self.dt_max_myr)?;
        positive("run.snapshot_interval_myr", self.snapshot_interval_myr)?;
        if whole_ratio(self.snapshot_interval_myr, self.dt_max_myr).is_none() {
            return err("run.snapshot_interval_myr must be a whole multiple of run.dt_max_myr");
        }
        Ok(())
    }
}

impl HydroConfig {
    pub(super) fn validate(&self) -> Result<(), ConfigError> {
        if self.gamma <= 1.0 {
            return err("hydro.gamma must exceed 1");
        }
        Ok(())
    }
}

impl BoxSection {
    pub(super) fn validate(&self) -> Result<(), ConfigError> {
        positive("box.side_kpc", self.side_kpc)
    }
}

impl GalaRPConfig {
    /// Checks `[galarp]` on its own. Whether `dt_myr` divides `run.dt_max_myr` is only
    /// checked by `Config::rk4_steps_per_dt_max`, as before.
    pub(super) fn validate(&self) -> Result<(), ConfigError> {
        positive("galarp.dt_myr", self.dt_myr)?;
        if let Some(s) = self.vertical_dispersion_km_s && !(s >= 0.0) {
            return err("galarp.vertical_dispersion_km_s must be >= 0");
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Times in code units
// ---------------------------------------------------------------------------

impl Config {
    pub fn t_end(&self) -> f64 {
        self.run.t_end_myr * MYR
    }

    pub fn dt_max(&self) -> f64 {
        self.run.dt_max_myr * MYR
    }

    /// Number of `dt_max` steps between snapshots.
    pub fn snapshot_every(&self) -> u64 {
        whole_ratio(self.run.snapshot_interval_myr, self.run.dt_max_myr).expect("validated")
    }

    /// RK4 steps per `run.dt_max_myr` (its log interval).
    pub fn rk4_steps_per_dt_max(&self) -> Result<u64, ConfigError> {
        match whole_ratio(self.run.dt_max_myr, self.galarp.dt_myr) {
            Some(n) => Ok(n),
            None => err("galarp.dt_myr must divide run.dt_max_myr"),
        }
    }
}
