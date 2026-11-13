//! `[ram_pressure]` and `[ram_pressure.shadow]`: the analytic wind on the clouds.

use galarp::analytic::{DensityAt, DynamicShadow, RamPressure};
use galarp::units::MSUN;
use serde::Deserialize;

use super::{err, positive, Config, ConfigError, ParamSpec};

// ---------------------------------------------------------------------------
// Schema
// ---------------------------------------------------------------------------

/// `[ram_pressure]`. Off unless `enabled = true`.
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct RamPressureConfig {
    pub enabled: bool,
    /// Cloud surface density Sigma: a number or a table in time.
    pub surface_density_msun_pc2: ParamSpec,
    /// Where rho is evaluated: "galaxy" (its center, as GalaRP) or "particle".
    pub density_at: String,
    pub shadow: ShadowConfig,
    /// Drag coefficient (efficiency of wind pushing on each cloud).
    pub drag_coefficient: f64,
}

/// `[ram_pressure.shadow]`: `DynamicShadow`. Off unless `enabled = true`.
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct ShadowConfig {
    pub enabled: bool,
    /// Side of the square columns along the wind.
    pub bin_width_kpc: f64,
    pub refill_length_kpc: f64,
    pub cutoff_radius_kpc: f64,
    pub cutoff_taper_kpc: f64,
}

// ---------------------------------------------------------------------------
// Defaults
// ---------------------------------------------------------------------------

impl Default for RamPressureConfig {
    fn default() -> Self {
        Self {
            enabled: false,
            surface_density_msun_pc2: ParamSpec::Value(10.0),
            density_at: "galaxy".into(),
            shadow: ShadowConfig::default(),
            drag_coefficient: 0.5,
        }
    }
}

impl Default for ShadowConfig {
    fn default() -> Self {
        Self { enabled: false,  
            bin_width_kpc: 2.0, 
            refill_length_kpc: f64::INFINITY, 
            cutoff_radius_kpc: f64::INFINITY, 
            cutoff_taper_kpc: 0.0 
        }
    }
}

// ---------------------------------------------------------------------------
// Validation (only when enabled; that it needs ICM gas is checked in `Config::validate`)
// ---------------------------------------------------------------------------

impl RamPressureConfig {
    pub(super) fn validate(&self) -> Result<(), ConfigError> {
        if !self.enabled {
            return Ok(());
        }
        let sigma = &self.surface_density_msun_pc2;
        sigma.validate("ram_pressure.surface_density_msun_pc2")?;
        let sigma_ok = match sigma {
            ParamSpec::Value(v) => *v > 0.0,
            ParamSpec::Table { value, .. } => value.iter().all(|&v| v > 0.0),
        };
        if !sigma_ok {
            return err("ram_pressure.surface_density_msun_pc2 must be positive");
        }
        if !["galaxy", "particle"].contains(&self.density_at.as_str()) {
            return err(format!("ram_pressure.density_at must be \"galaxy\" or \"particle\", got \"{}\"", self.density_at));
        }
        if self.shadow.enabled {
            positive("ram_pressure.shadow.bin_width_kpc", self.shadow.bin_width_kpc)?;
        }
        
        positive("ram_pressure.drag_coefficient", self.drag_coefficient)?;
        
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Conversion to code units / galarp types
// ---------------------------------------------------------------------------

impl Config {
    /// The ram-pressure model, or `None` if `[ram_pressure]` is off.
    pub fn ram_pressure(&self) -> Result<Option<RamPressure>, ConfigError> {
        let rp = &self.ram_pressure;
        if !rp.enabled {
            return Ok(None);
        }
        // Msun / pc^2 to 1e10 Msun / kpc^2.
        let sigma = rp.surface_density_msun_pc2.to_param(1e6 * MSUN);
        let drag_coeff = rp.drag_coefficient;
        let shadow = rp.shadow.enabled.then_some(DynamicShadow {
            bin_width: rp.shadow.bin_width_kpc,
            refill_length: rp.shadow.refill_length_kpc,
            cutoff_radius: rp.shadow.cutoff_radius_kpc,
            cutoff_taper: rp.shadow.cutoff_taper_kpc,
        });

        let density_at = if rp.density_at == "particle" { DensityAt::Particle } else { DensityAt::Galaxy };
        Ok(Some(RamPressure { 
            icm: self.cluster_gas()?, 
            surface_density: sigma, 
            density_at, 
            drag_coeff, 
            shadow 
        }))
    }
}
