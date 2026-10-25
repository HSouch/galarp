
pub mod beta_profile;
pub mod hse;
pub mod tabulated;

pub use beta_profile::BetaModel;
pub use tabulated::RadialTable;

use crate::potentials::Potential;


pub enum ClusterGas {
    Uniform{rho: f64, u: f64},
    Radial(RadialTable),
}

impl ClusterGas {

    pub fn beta_model_hse(
        model: &BetaModel, 
        potential: &dyn Potential,
        gamma: f64,
        r_min: f64,
        r_max: f64,
        n: usize,
    ) -> Self {
        ClusterGas::Radial(
            hse::hydrostatic_table(
                |r| model.density(r), 
                potential, 
                gamma, 
                r_min, 
                r_max, 
                n))
    }


    pub fn density(&self, r: f64) -> f64 {
        match self {
            ClusterGas::Uniform { rho, .. } => *rho,
            ClusterGas::Radial(t) => t.density(r),
        }
    }

    pub fn specific_energy(&self, r: f64) -> f64 {
        match self {
            ClusterGas::Uniform { u, .. } => *u,
            ClusterGas::Radial(t) => t.specific_energy(r),
        }
    }

    /// Density at a cluster-frame position.
    pub fn density_at(&self, x: &[f64; 3]) -> f64 {
        self.density((x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt())
    }

    /// Specific internal energy at a cluster-frame position.
    pub fn specific_energy_at(&self, x: &[f64; 3]) -> f64 {
        self.specific_energy((x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt())
    }
}