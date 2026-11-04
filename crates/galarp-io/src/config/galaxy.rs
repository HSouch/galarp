//! `[[galaxy]]` sections: each galaxy's potential, orientation, gas disk and orbit.

use galarp::ic::{DiskOptions, GasDiskParams};
use galarp::orientation::Orientation;
use galarp::potentials::Composite;
use galarp::units::MSUN;
use serde::Deserialize;

use super::{composite, err, positive, ConfigError, OrbitConfig, PotentialTerm};

// ---------------------------------------------------------------------------
// Schema
// ---------------------------------------------------------------------------

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GalaxyConfig {
    /// Optional label, used in the log and printed summaries.
    pub name: Option<String>,
    /// Static galaxy potential terms (dark matter + stars).
    pub potential: Vec<PotentialTerm>,
    /// Orientation of the galaxy (its disk and potential): the disk's spin axis in the
    /// cluster frame, either as a vector or as an inclination from +z tipped toward an
    /// azimuth in the x-y plane (from +x toward +y). Default: [0, 0, 1], the x-y plane.
    pub disk_normal: Option<[f64; 3]>,
    pub inclination_deg: Option<f64>,
    pub azimuth_deg: Option<f64>,
    /// The gas disk the clouds are sampled from. Without one the galaxy only gravitates.
    pub disk: Option<DiskConfig>,
    pub orbit: OrbitConfig,
}

/// `[galaxy.disk]`
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DiskConfig {
    pub mass_msun: f64,
    pub scale_length_kpc: f64,
    pub r_max_kpc: f64,
    /// Effective isothermal sound speed (thermal plus turbulent support). Sets the disk's
    /// vertical profile and the clouds' default vertical dispersion.
    pub sound_speed_km_s: f64,
    /// Number of clouds; each carries mass_msun / n_particles.
    pub n_particles: usize,
    /// Include the gas disk's own gravity in its vertical structure (and so in where the
    /// clouds start).
    #[serde(default)]
    pub include_self_gravity: bool,
    /// Gaussian position noise, as a fraction of the local cloud spacing (0 = none).
    #[serde(default)]
    pub noise_fraction: f64,
    /// Seed for that noise.
    #[serde(default = "d_noise_seed")]
    pub noise_seed: u64,
}

// ---------------------------------------------------------------------------
// Defaults
// ---------------------------------------------------------------------------

fn d_noise_seed() -> u64 {
    1
}

// ---------------------------------------------------------------------------
// Validation (this galaxy only; the caller adds the "galaxy k (name):" prefix)
// ---------------------------------------------------------------------------

impl GalaxyConfig {
    pub(super) fn validate(&self) -> Result<(), ConfigError> {
        if let Some(d) = &self.disk {
            d.validate()?;
        }
        for t in &self.potential {
            t.validate()?;
        }
        self.validate_orientation()
    }

    fn validate_orientation(&self) -> Result<(), ConfigError> {
        if self.disk_normal.is_some() && (self.inclination_deg.is_some() || self.azimuth_deg.is_some()) {
            return err("give disk_normal or inclination_deg/azimuth_deg, not both");
        }
        if self.azimuth_deg.is_some() && self.inclination_deg.is_none() {
            return err("azimuth_deg needs inclination_deg");
        }
        if let Some(n) = self.disk_normal {
            if !(n.iter().all(|c| c.is_finite()) && n.iter().any(|&c| c != 0.0)) {
                return err("disk_normal must be a nonzero vector");
            }
        }
        if let Some(i) = self.inclination_deg {
            if !(0.0..=180.0).contains(&i) {
                return err("inclination_deg must be between 0 and 180");
            }
        }
        Ok(())
    }
}

impl DiskConfig {
    fn validate(&self) -> Result<(), ConfigError> {
        positive("disk.mass_msun", self.mass_msun)?;
        positive("disk.scale_length_kpc", self.scale_length_kpc)?;
        positive("disk.r_max_kpc", self.r_max_kpc)?;
        positive("disk.sound_speed_km_s", self.sound_speed_km_s)?;
        if self.n_particles == 0 {
            return err("disk.n_particles must be at least 1");
        }
        if !(self.noise_fraction >= 0.0) {
            return err("disk.noise_fraction must be >= 0");
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Conversion to code units / galarp types
// ---------------------------------------------------------------------------

impl GalaxyConfig {
    /// Its `name`, or "galaxy k" (`k` is its index in the config).
    pub fn label(&self, k: usize) -> String {
        self.name.clone().unwrap_or_else(|| format!("galaxy {k}"))
    }

    /// Orientation of the galaxy (its disk and potential).
    pub fn orientation(&self) -> Orientation {
        match (self.disk_normal, self.inclination_deg) {
            (Some(n), _) => Orientation::from_normal(n),
            (None, Some(i)) => Orientation::from_angles(i, self.azimuth_deg.unwrap_or(0.0)),
            (None, None) => Orientation::IDENTITY,
        }
    }

    /// The galaxy's potential, in code units.
    pub fn build_potential(&self) -> Composite {
        composite(&self.potential)
    }

    /// Initial-condition options for the disk (self-gravity, noise); defaults without one.
    pub fn disk_options(&self) -> DiskOptions {
        self.disk.as_ref().map(DiskConfig::options).unwrap_or_default()
    }
}

impl DiskConfig {
    /// Disk parameters in code units.
    pub fn params(&self) -> GasDiskParams {
        GasDiskParams {
            mass: self.mass_msun * MSUN,
            scale_length: self.scale_length_kpc,
            r_max: self.r_max_kpc,
            sound_speed: self.sound_speed_km_s,
        }
    }

    pub fn options(&self) -> DiskOptions {
        DiskOptions {
            self_gravity: self.include_self_gravity,
            noise_fraction: self.noise_fraction,
            noise_seed: self.noise_seed,
        }
    }
}
