//! `[galaxy.orbit]` (how each galaxy moves) and `[orbits]` (how integrated orbits are
//! coupled), plus building every galaxy's orbit table.

use std::path::{Path, PathBuf};

use galarp::orbit::{self, Body, Motion, OrbitTable};
use galarp::potentials::{Composite, Potential};
use galarp::units::MYR;
use serde::Deserialize;

use super::{err, Config, ConfigError};

// ---------------------------------------------------------------------------
// Schema
// ---------------------------------------------------------------------------

/// `[galaxy.orbit]`
#[derive(Clone, Debug, Deserialize)]
#[serde(tag = "mode", rename_all = "snake_case", deny_unknown_fields)]
pub enum OrbitConfig {
    /// Straight line from `position_kpc` (cluster frame) at t = 0.
    ConstantVelocity { position_kpc: [f64; 3], velocity_km_s: [f64; 3] },
    /// Orbit table CSV as written by `OrbitTable::write_csv` (code units: kpc, km/s,
    /// kpc/(km/s)).
    Table { file: PathBuf },
    /// Integrate the galaxy as a point mass in the cluster potential from t = 0.
    Integrate {
        position_kpc: [f64; 3],
        velocity_km_s: [f64; 3],
        #[serde(default = "d_orbit_dt")]
        dt_myr: f64,
        #[serde(default = "d_orbit_out")]
        output_dt_myr: f64,
    },
}

/// `[orbits]`: how galaxies with `orbit.mode = "integrate"` move.
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct OrbitsConfig {
    /// Each integrated galaxy is also pulled by every other galaxy's potential (see
    /// `galarp::orbit::mutual`). Off: each feels only the cluster potential.
    pub mutual_gravity: bool,
}

// ---------------------------------------------------------------------------
// Defaults
// ---------------------------------------------------------------------------

fn d_orbit_dt() -> f64 {
    0.1
}

fn d_orbit_out() -> f64 {
    1.0
}

impl Default for OrbitsConfig {
    fn default() -> Self {
        Self { mutual_gravity: true }
    }
}

// ---------------------------------------------------------------------------
// Conversion to code units / galarp types
// ---------------------------------------------------------------------------

impl OrbitConfig {
    pub fn is_integrated(&self) -> bool {
        matches!(self, OrbitConfig::Integrate { .. })
    }

    /// `(dt, output_dt)` of an `integrate` orbit, in code units.
    fn integration_steps(&self) -> Option<(f64, f64)> {
        match self {
            OrbitConfig::Integrate { dt_myr, output_dt_myr, .. } => Some((dt_myr * MYR, output_dt_myr * MYR)),
            _ => None,
        }
    }

    /// The galaxy's starting `Motion`: a prescribed table that covers [0, t_end] (a
    /// constant-velocity line is extended to `t_cover`), or the initial state to integrate.
    /// `label` names the galaxy in error messages.
    fn motion(&self, label: &str, t_end: f64, t_cover: f64) -> Result<Motion, ConfigError> {
        match self {
            OrbitConfig::ConstantVelocity { position_kpc: x, velocity_km_s: v } => {
                let at = |t: f64| [x[0] + v[0] * t, x[1] + v[1] * t, x[2] + v[2] * t];
                let t1 = (1.01 * t_end).max(t_cover);
                let table = OrbitTable::new(vec![0.0, t1], vec![at(0.0), at(t1)], vec![*v, *v])
                    .map_err(|e| ConfigError(e.to_string()))?;
                Ok(Motion::Prescribed(table))
            }
            OrbitConfig::Table { file } => {
                let table = read_orbit_csv(file)?;
                if table.t_start() > 0.0 || table.t_end() < t_end {
                    return err(format!(
                        "{label}: orbit covers t = {:.1}-{:.1} Myr but the run needs 0-{:.1} Myr",
                        table.t_start() / MYR,
                        table.t_end() / MYR,
                        t_end / MYR
                    ));
                }
                Ok(Motion::Prescribed(table))
            }
            OrbitConfig::Integrate { position_kpc, velocity_km_s, .. } => {
                Ok(Motion::Integrate { x0: *position_kpc, v0: *velocity_km_s })
            }
        }
    }
}

fn read_orbit_csv(file: &Path) -> Result<OrbitTable, ConfigError> {
    let in_file = |e: &dyn std::fmt::Display| ConfigError(format!("{}: {e}", file.display()));
    let f = std::fs::File::open(file).map_err(|e| in_file(&e))?;
    OrbitTable::read_csv(std::io::BufReader::new(f)).map_err(|e| in_file(&e))
}

impl Config {
    /// Orbit of galaxy `k` in code units, covering at least [0, t_end]. Computes all of
    /// them (see `orbit_tables`); call that once when you need several.
    pub fn orbit_table(&self, k: usize) -> Result<OrbitTable, ConfigError> {
        Ok(self.orbit_tables()?.swap_remove(k))
    }

    /// Every galaxy's orbit in code units, each covering at least [0, t_end].
    ///
    /// Galaxies with `constant_velocity` or `table` orbits follow them. Galaxies with
    /// `integrate` orbits are integrated together (`orbit::integrate_bodies`), pulled by the
    /// cluster potential and, with `orbits.mutual_gravity`, by every other galaxy's
    /// potential, including the prescribed ones. They share the smallest `dt_myr` and
    /// `output_dt_myr` among them.
    pub fn orbit_tables(&self) -> Result<Vec<OrbitTable>, ConfigError> {
        let t_end = self.t_end();

        // 1. Shared step sizes of the integrated orbits. They run one output step past
        //    t_end so the table still covers it after rounding.
        let steps: Vec<(f64, f64)> = self.galaxy.iter().filter_map(|g| g.orbit.integration_steps()).collect();
        let dt = steps.iter().map(|s| s.0).fold(f64::INFINITY, f64::min);
        let output_dt = steps.iter().map(|s| s.1).fold(f64::INFINITY, f64::min);
        let t_integrate = if steps.is_empty() { t_end } else { t_end + output_dt };

        // 2. Each galaxy's starting motion.
        let motions: Vec<Motion> = self
            .galaxy
            .iter()
            .enumerate()
            .map(|(k, g)| g.orbit.motion(&g.label(k), t_end, t_integrate))
            .collect::<Result<_, _>>()?;

        // 3. Nothing to integrate: every orbit is already a table.
        if steps.is_empty() {
            return Ok(motions
                .into_iter()
                .map(|m| match m {
                    Motion::Prescribed(t) => t,
                    Motion::Integrate { .. } => unreachable!(),
                })
                .collect());
        }

        // 4. Integrate. Prescribed tables must reach the end of the integration; stop it
        //    early if one does not (it still covers t_end, checked in `motion`).
        let t1 = motions
            .iter()
            .filter_map(|m| match m {
                Motion::Prescribed(t) => Some(t.t_end()),
                Motion::Integrate { .. } => None,
            })
            .fold(t_integrate, f64::min);
        let potentials: Vec<Composite> = self.galaxy.iter().map(|g| g.build_potential()).collect();
        let bodies: Vec<Body> = motions
            .into_iter()
            .zip(&self.galaxy)
            .zip(&potentials)
            .map(|((motion, g), potential)| Body { potential, orientation: g.orientation(), motion })
            .collect();
        let cluster = self.cluster_potential();
        let cluster_ref = cluster.as_ref().map(|c| c as &dyn Potential);
        Ok(orbit::integrate_bodies(cluster_ref, &bodies, self.orbits.mutual_gravity, 0.0, t1, dt, output_dt.min(t1)))
    }
}
