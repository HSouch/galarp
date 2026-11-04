//! `[[galaxy.potential]]` and `[[cluster.potential]]` terms.

use galarp::potentials::{Burkert, Composite, Hernquist, MiyamotoNagai, Nfw, PointMass, Potential};
use galarp::units::MSUN;
use serde::Deserialize;

use super::{err, ConfigError, ParamSpec};

// ---------------------------------------------------------------------------
// Schema
// ---------------------------------------------------------------------------

/// One potential term. `type` picks the profile; every `ParamSpec` field can be a number
/// or a table in time.
///
/// To add a profile: add a variant here, then one arm each in `kind`, `params` and `build`.
#[derive(Clone, Debug, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case", deny_unknown_fields)]
pub enum PotentialTerm {
    PointMass { mass_msun: ParamSpec },
    Hernquist { mass_msun: ParamSpec, a_kpc: ParamSpec },
    MiyamotoNagai { mass_msun: ParamSpec, a_kpc: ParamSpec, b_kpc: ParamSpec },
    /// NFW with Gala's scale mass m = 4 pi rho_0 r_s^3.
    Nfw { m_msun: ParamSpec, r_s_kpc: ParamSpec },
    /// NFW from virial mass, virial radius and concentration (constant in time).
    NfwVirial { m_vir_msun: f64, r_vir_kpc: f64, concentration: f64 },
    Burkert { rho0_msun_kpc3: ParamSpec, r0_kpc: ParamSpec },
}

impl PotentialTerm {
    /// The term's `type` as written in the config.
    pub fn kind(&self) -> &'static str {
        match self {
            PotentialTerm::PointMass { .. } => "point_mass",
            PotentialTerm::Hernquist { .. } => "hernquist",
            PotentialTerm::MiyamotoNagai { .. } => "miyamoto_nagai",
            PotentialTerm::Nfw { .. } => "nfw",
            PotentialTerm::NfwVirial { .. } => "nfw_virial",
            PotentialTerm::Burkert { .. } => "burkert",
        }
    }

    /// The term's `ParamSpec` fields with their config keys. (`nfw_virial` has none: its
    /// fields are plain numbers.)
    fn params(&self) -> Vec<(&'static str, &ParamSpec)> {
        match self {
            PotentialTerm::PointMass { mass_msun } => vec![("mass_msun", mass_msun)],
            PotentialTerm::Hernquist { mass_msun, a_kpc } => vec![("mass_msun", mass_msun), ("a_kpc", a_kpc)],
            PotentialTerm::MiyamotoNagai { mass_msun, a_kpc, b_kpc } => {
                vec![("mass_msun", mass_msun), ("a_kpc", a_kpc), ("b_kpc", b_kpc)]
            }
            PotentialTerm::Nfw { m_msun, r_s_kpc } => vec![("m_msun", m_msun), ("r_s_kpc", r_s_kpc)],
            PotentialTerm::NfwVirial { .. } => vec![],
            PotentialTerm::Burkert { rho0_msun_kpc3, r0_kpc } => vec![("rho0_msun_kpc3", rho0_msun_kpc3), ("r0_kpc", r0_kpc)],
        }
    }
}

// ---------------------------------------------------------------------------
// Validation
// ---------------------------------------------------------------------------

impl PotentialTerm {
    pub(super) fn validate(&self) -> Result<(), ConfigError> {
        for (key, p) in self.params() {
            p.validate(&format!("{}.{key}", self.kind()))?;
        }
        if let PotentialTerm::NfwVirial { m_vir_msun, r_vir_kpc, concentration } = self {
            if !(*m_vir_msun > 0.0 && *r_vir_kpc > 0.0 && *concentration > 0.0) {
                return err("nfw_virial: m_vir_msun, r_vir_kpc and concentration must be positive");
            }
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Conversion to code units / galarp types
// ---------------------------------------------------------------------------

impl PotentialTerm {
    /// The term as a galarp potential, in code units (masses in 1e10 Msun, lengths in kpc).
    pub fn build(&self) -> Box<dyn Potential> {
        match self {
            PotentialTerm::PointMass { mass_msun } => Box::new(PointMass { mass: mass_msun.to_param(MSUN) }),
            PotentialTerm::Hernquist { mass_msun, a_kpc } => {
                Box::new(Hernquist { mass: mass_msun.to_param(MSUN), a: a_kpc.to_param(1.0) })
            }
            PotentialTerm::MiyamotoNagai { mass_msun, a_kpc, b_kpc } => Box::new(MiyamotoNagai {
                mass: mass_msun.to_param(MSUN),
                a: a_kpc.to_param(1.0),
                b: b_kpc.to_param(1.0),
            }),
            PotentialTerm::Nfw { m_msun, r_s_kpc } => {
                Box::new(Nfw { m: m_msun.to_param(MSUN), r_s: r_s_kpc.to_param(1.0) })
            }
            PotentialTerm::NfwVirial { m_vir_msun, r_vir_kpc, concentration } => {
                Box::new(Nfw::from_virial(m_vir_msun * MSUN, *r_vir_kpc, *concentration))
            }
            PotentialTerm::Burkert { rho0_msun_kpc3, r0_kpc } => {
                Box::new(Burkert { rho_0: rho0_msun_kpc3.to_param(MSUN), r_0: r0_kpc.to_param(1.0) })
            }
        }
    }
}

/// The sum of `terms`.
pub(super) fn composite(terms: &[PotentialTerm]) -> Composite {
    Composite { terms: terms.iter().map(|t| t.build()).collect() }
}
