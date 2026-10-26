
//! TOML Configuration files for GalaRP
//! 
//! Reading, writing, and validation


use std::fmt;
use std::path::{Path, PathBuf};

use bigwig::cluster::{hse, BetaModel, ClusterGas, RadialTable};
use bigwig::orbit::{self, Body, Motion, OrbitTable};
use bigwig::potentials::{Burkert, Composite, Hernquist, MiyamotoNagai, Nfw, Param, PointMass, Potential};
use bigwig::units::{self, UNIT_DENSITY_CGS, MSUN, MYR};
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



#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RunConfig {
    /// Output directory.. 
    /// Relative paths are from the directory the program is run from.
    pub output_dir: PathBuf,
    pub t_end_myr: f64,
    /// Log interval`.
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


#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct HydroConfig {
    pub gamma: f64,
}

impl Default for HydroConfig {
    fn default() -> Self {
        Self { gamma: units::GAMMA }
    }
}


#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BoxSection {
    pub side_kpc: f64,
}

#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ClusterConfig {
    pub gas: ClusterGasConfig,
    /// Cluster potential terms (empty: no cluster gravity, as for a uniform wind).
    #[serde(default)]
    pub potential: Vec<PotentialTerm>,
}


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

fn d_table_rmin() -> f64 {
    1.0
}
fn d_table_rmax() -> f64 {
    1e5
}
fn d_table_n() -> usize {
    4000
}
fn d_orbit_dt() -> f64 {
    0.1
}
fn d_orbit_out() -> f64 {
    1.0
}


/// A potential parameter: a number, or a table in time.
#[derive(Clone, Debug, Deserialize)]
#[serde(untagged)]
pub enum ParamSpec {
    Value(f64),
    Table { t_myr: Vec<f64>, value: Vec<f64> },
}

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

fn d_noise_seed() -> u64 {
    1
}


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


/// Region counted as "still in the disk" for the log's `in_disk_fraction` column: a
/// cylinder around the galaxy center, aligned with the disk.
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct DiagnosticsConfig {
    pub disk_region_radius_kpc: f64,
    pub disk_region_half_height_kpc: f64,
}

impl Default for DiagnosticsConfig {
    fn default() -> Self {
        Self { disk_region_radius_kpc: 15.0, disk_region_half_height_kpc: 2.0 }
    }
}


/// How galaxies with `orbit.mode = "integrate"` move.
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct OrbitsConfig {
    /// Each integrated galaxy is also pulled by every other galaxy's potential (see
    /// `bigwig::orbit::mutual`). Off: each feels only the cluster potential.
    pub mutual_gravity: bool,
}

impl Default for OrbitsConfig {
    fn default() -> Self {
        Self { mutual_gravity: true }
    }
}


#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct RamPressureConfig {
    pub enabled: bool,
    /// Cloud surface density Sigma: a number or a table in time.
    pub surface_density_msun_pc2: ParamSpec,
    /// Where rho is evaluated: "galaxy" (its center, as GalaRP) or "particle".
    pub density_at: String,
    pub shadow: ShadowConfig,
}

impl Default for RamPressureConfig {
    fn default() -> Self {
        Self { enabled: false, 
            surface_density_msun_pc2: ParamSpec::Value(10.0), density_at: "galaxy".into(), 
            shadow: ShadowConfig::default() }
    }
}


/// `DynamicShadow` Off unless `enabled = true`.
#[derive(Clone, Debug, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct ShadowConfig {
    pub enabled: bool,
    /// C_D in the optical depth tau = C_D A / bin_width^2 of each particle.
    pub drag_coefficient: f64,
    /// Side of the square columns along the wind.
    pub bin_width_kpc: f64,
}

impl Default for ShadowConfig {
    fn default() -> Self {
        Self { enabled: false, drag_coefficient: 0.5, bin_width_kpc: 2.0 }
    }
}


/// GalaRP-specific Config Flags
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

impl Default for GalaRPConfig {
    fn default() -> Self {
        Self { dt_myr: 0.5, vertical_dispersion_km_s: None, seed: 1 }
    }
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
        let positive = |name: &str, v: f64| if v > 0.0 && v.is_finite() { Ok(()) } else { err(format!("{name} must be positive, got {v}")) };
        positive("run.t_end_myr", self.run.t_end_myr)?;
        positive("run.dt_max_myr", self.run.dt_max_myr)?;
        positive("run.snapshot_interval_myr", self.run.snapshot_interval_myr)?;
        positive("box.side_kpc", self.box_.side_kpc)?;
        let ratio = self.run.snapshot_interval_myr / self.run.dt_max_myr;
        if (ratio - ratio.round()).abs() > 1e-9 || ratio.round() < 1.0 {
            return err("run.snapshot_interval_myr must be a whole multiple of run.dt_max_myr");
        }
        if self.hydro.gamma <= 1.0 {
            return err("hydro.gamma must exceed 1");
        }
        if self.galaxy.is_empty() {
            return err("at least one [[galaxy]] is required; the first is the parent");
        }
        if self.galaxy.len() > u16::MAX as usize {
            return err("too many galaxies");
        }
        for (k, g) in self.galaxy.iter().enumerate() {
            let at = |key: &str| format!("galaxy {} ({}): {key}", k, self.galaxy_name(k));
            if let Some(d) = &g.disk {
                positive(&at("disk.mass_msun"), d.mass_msun)?;
                positive(&at("disk.scale_length_kpc"), d.scale_length_kpc)?;
                positive(&at("disk.r_max_kpc"), d.r_max_kpc)?;
                positive(&at("disk.sound_speed_km_s"), d.sound_speed_km_s)?;
                if d.n_particles == 0 {
                    return err(at("disk.n_particles must be at least 1"));
                }
                if !(d.noise_fraction >= 0.0) {
                    return err(at("disk.noise_fraction must be >= 0"));
                }
            }
            for t in &g.potential {
                t.validate().map_err(|e| ConfigError(at(&e.0)))?;
            }
            if g.disk_normal.is_some() && (g.inclination_deg.is_some() || g.azimuth_deg.is_some()) {
                return err(at("give disk_normal or inclination_deg/azimuth_deg, not both"));
            }
            if g.azimuth_deg.is_some() && g.inclination_deg.is_none() {
                return err(at("azimuth_deg needs inclination_deg"));
            }
            if let Some(n) = g.disk_normal {
                if !(n.iter().all(|c| c.is_finite()) && n.iter().any(|&c| c != 0.0)) {
                    return err(at("disk_normal must be a nonzero vector"));
                }
            }
            if let Some(i) = g.inclination_deg {
                if !(0.0..=180.0).contains(&i) {
                    return err(at("inclination_deg must be between 0 and 180"));
                }
            }
        }
        // An integrated orbit needs something to pull on it: the cluster, or other galaxies.
        let pulled_by_galaxies = self.orbits.mutual_gravity && self.galaxy.len() > 1;
        let needs_cluster_potential = matches!(self.cluster.gas, ClusterGasConfig::Beta { .. } | ClusterGasConfig::Table { hydrostatic: true, .. })
            || (!pulled_by_galaxies && self.galaxy.iter().any(|g| matches!(g.orbit, OrbitConfig::Integrate { .. })));
        if needs_cluster_potential && self.cluster.potential.is_empty() {
            return err("this cluster gas model or orbit mode needs [[cluster.potential]] terms");
        }
        for t in &self.cluster.potential {
            t.validate()?;
        }
        let rp = &self.ram_pressure;
        if rp.enabled {
            rp.surface_density_msun_pc2.validate("ram_pressure.surface_density_msun_pc2")?;
            let sigma_ok = match &rp.surface_density_msun_pc2 {
                ParamSpec::Value(v) => *v > 0.0,
                ParamSpec::Table { value, .. } => value.iter().all(|&v| v > 0.0),
            };
            if !sigma_ok {
                return err("ram_pressure.surface_density_msun_pc2 must be positive");
            }
            if !["galaxy", "particle"].contains(&rp.density_at.as_str()) {
                return err(format!("ram_pressure.density_at must be \"galaxy\" or \"particle\", got \"{}\"", rp.density_at));
            }
            if matches!(self.cluster.gas, ClusterGasConfig::None) {
                return err("ram_pressure needs ICM gas: [cluster.gas] model = \"none\" has no density");
            }
            if rp.shadow.enabled {
                positive("ram_pressure.shadow.drag_coefficient", rp.shadow.drag_coefficient)?;
                positive("ram_pressure.shadow.bin_width_kpc", rp.shadow.bin_width_kpc)?;
            }
        }
        let r = &self.galarp;
        positive("galarp.dt_myr", r.dt_myr)?;
        if let Some(s) = r.vertical_dispersion_km_s {
            if !(s >= 0.0) {
                return err("galarp.vertical_dispersion_km_s must be >= 0");
            }
        }
        Ok(())
    }
}


impl fmt::Display for ParamSpec {
    /// The value, or "table" for a table in time.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParamSpec::Value(v) => write!(f, "{v}"),
            ParamSpec::Table { .. } => f.write_str("(table in time)"),
        }
    }
}

impl ParamSpec {
    fn validate(&self, name: &str) -> Result<(), ConfigError> {
        if let ParamSpec::Table { t_myr, value } = self {
            if t_myr.is_empty() || t_myr.len() != value.len() {
                return err(format!("{name}: t_myr and value must be non-empty and the same length"));
            }
            if !t_myr.windows(2).all(|w| w[0] < w[1]) {
                return err(format!("{name}: t_myr must increase"));
            }
        }
        Ok(())
    }

    /// Converts to a code-unit `Param`, multiplying values by `scale`.
    pub fn to_param(&self, scale: f64) -> Param {
        match self {
            ParamSpec::Value(v) => Param::Constant(v * scale),
            ParamSpec::Table { t_myr, value } => {
                Param::table(t_myr.iter().map(|t| t * MYR).collect(), value.iter().map(|v| v * scale).collect())
            }
        }
    }
}


impl PotentialTerm {
    fn validate(&self) -> Result<(), ConfigError> {
        match self {
            PotentialTerm::PointMass { mass_msun } => mass_msun.validate("point_mass.mass_msun"),
            PotentialTerm::Hernquist { mass_msun, a_kpc } => {
                mass_msun.validate("hernquist.mass_msun")?;
                a_kpc.validate("hernquist.a_kpc")
            }
            PotentialTerm::MiyamotoNagai { mass_msun, a_kpc, b_kpc } => {
                mass_msun.validate("miyamoto_nagai.mass_msun")?;
                a_kpc.validate("miyamoto_nagai.a_kpc")?;
                b_kpc.validate("miyamoto_nagai.b_kpc")
            }
            PotentialTerm::Nfw { m_msun, r_s_kpc } => {
                m_msun.validate("nfw.m_msun")?;
                r_s_kpc.validate("nfw.r_s_kpc")
            }
            PotentialTerm::NfwVirial { m_vir_msun, r_vir_kpc, concentration } => {
                if *m_vir_msun > 0.0 && *r_vir_kpc > 0.0 && *concentration > 0.0 {
                    Ok(())
                } else {
                    err("nfw_virial: m_vir_msun, r_vir_kpc and concentration must be positive")
                }
            }
            PotentialTerm::Burkert { rho0_msun_kpc3, r0_kpc } => {
                rho0_msun_kpc3.validate("burkert.rho0_msun_kpc3")?;
                r0_kpc.validate("burkert.r0_kpc")
            }
        }
    }

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
            PotentialTerm::Nfw { m_msun, r_s_kpc } => Box::new(Nfw { m: m_msun.to_param(MSUN), r_s: r_s_kpc.to_param(1.0) }),
            PotentialTerm::NfwVirial { m_vir_msun, r_vir_kpc, concentration } => {
                Box::new(Nfw::from_virial(m_vir_msun * MSUN, *r_vir_kpc, *concentration))
            }
            PotentialTerm::Burkert { rho0_msun_kpc3, r0_kpc } => {
                Box::new(Burkert { rho_0: rho0_msun_kpc3.to_param(MSUN), r_0: r0_kpc.to_param(1.0) })
            }
        }
    }
}

fn composite(terms: &[PotentialTerm]) -> Composite {
    Composite { terms: terms.iter().map(|t| t.build()).collect() }
}


impl Config {
    pub fn ram_pressure(&self) -> Result<Option<bigwig::galarp::RamPressure>, ConfigError> {
        use bigwig::galarp::{DensityAt, DynamicShadow, RamPressure};
        let rp = &self.ram_pressure;
        if !rp.enabled {
            return Ok(None);
        }
        // Msun / pc^2 to 1e10 Msun / kpc^2.
        let sigma = rp.surface_density_msun_pc2.to_param(1e6 * MSUN);
        let shadow = rp.shadow.enabled.then_some(DynamicShadow { drag_coeff: rp.shadow.drag_coefficient, bin_width: rp.shadow.bin_width_kpc });
        let density_at = if rp.density_at == "particle" { DensityAt::Particle } else { DensityAt::Galaxy };
        Ok(Some(RamPressure { icm: self.cluster_gas()?, surface_density: sigma, density_at, shadow }))
    }

    /// RK4 steps per  `run.dt_max_myr` (its log interval).
    pub fn rk4_steps_per_dt_max(&self) -> Result<u64, ConfigError> {
        let ratio = self.run.dt_max_myr / self.galarp.dt_myr;
        if (ratio - ratio.round()).abs() > 1e-9 || ratio.round() < 1.0 {
            return err("galarp.dt_myr must divide run.dt_max_myr");
        }
        Ok(ratio.round() as u64)
    }

    pub fn t_end(&self) -> f64 {
        self.run.t_end_myr * MYR
    }

    pub fn dt_max(&self) -> f64 {
        self.run.dt_max_myr * MYR
    }

    /// Number of `dt_max` steps between snapshots.
    pub fn snapshot_every(&self) -> u64 {
        (self.run.snapshot_interval_myr / self.run.dt_max_myr).round() as u64
    }

    /// Label of galaxy `k`: its `name`, or "galaxy k".
    pub fn galaxy_name(&self, k: usize) -> String {
        self.galaxy[k].name.clone().unwrap_or_else(|| format!("galaxy {k}"))
    }

    /// Initial-condition options for galaxy `k`'s disk (self-gravity, noise).
    pub fn disk_options(&self, k: usize) -> bigwig::ic::DiskOptions {
        let Some(d) = self.galaxy[k].disk.as_ref() else { return bigwig::ic::DiskOptions::default() };
        bigwig::ic::DiskOptions { self_gravity: d.include_self_gravity, noise_fraction: d.noise_fraction, noise_seed: d.noise_seed }
    }

    /// Orientation of galaxy `k` (its disk and potential).
    pub fn galaxy_orientation(&self, k: usize) -> bigwig::orientation::Orientation {
        use bigwig::orientation::Orientation;
        let g = &self.galaxy[k];
        match (g.disk_normal, g.inclination_deg) {
            (Some(n), _) => Orientation::from_normal(n),
            (None, Some(i)) => Orientation::from_angles(i, g.azimuth_deg.unwrap_or(0.0)),
            (None, None) => Orientation::IDENTITY,
        }
    }

    /// Disk parameters (code units) and cloud count of galaxy `k`, if it has a disk.
    pub fn disk(&self, k: usize) -> Option<(bigwig::ic::GasDiskParams, usize)> {
        let d = self.galaxy[k].disk.as_ref()?;
        let params = bigwig::ic::GasDiskParams {
            mass: d.mass_msun * MSUN,
            scale_length: d.scale_length_kpc,
            r_max: d.r_max_kpc,
            sound_speed: d.sound_speed_km_s,
        };
        Some((params, d.n_particles))
    }

    /// Cluster potential, or `None` if no terms are given.
    pub fn cluster_potential(&self) -> Option<Composite> {
        (!self.cluster.potential.is_empty()).then(|| composite(&self.cluster.potential))
    }

    /// Potential of galaxy `k` (0 = parent).
    pub fn galaxy_potential(&self, k: usize) -> Composite {
        composite(&self.galaxy[k].potential)
    }

    pub fn cluster_gas(&self) -> Result<ClusterGas, ConfigError> {
        let gamma = self.hydro.gamma;
        match &self.cluster.gas {
            ClusterGasConfig::None => Ok(ClusterGas::Uniform { rho: 0.0, u: 0.0 }),
            ClusterGasConfig::Uniform { density_g_cm3, temperature_k, mean_molecular_weight } => Ok(ClusterGas::Uniform {
                rho: density_g_cm3 / UNIT_DENSITY_CGS,
                u: units::u_from_temperature(*temperature_k, *mean_molecular_weight, gamma),
            }),
            ClusterGasConfig::Beta { rho0_g_cm3, r_c_kpc, beta, r_min_kpc, r_max_kpc, n_table } => {
                let model = BetaModel { rho_0: rho0_g_cm3 / UNIT_DENSITY_CGS, r_c: *r_c_kpc, beta: *beta };
                let pot = self.cluster_potential().expect("validated");
                Ok(ClusterGas::beta_model_hse(&model, &pot, gamma, *r_min_kpc, *r_max_kpc, *n_table))
            }
            ClusterGasConfig::Table { file, mean_molecular_weight, hydrostatic } => {
                let (r, rho, t) = read_profile(file)?;
                let rho: Vec<f64> = rho.iter().map(|v| v / UNIT_DENSITY_CGS).collect();
                if *hydrostatic {
                    let pot = self.cluster_potential().expect("validated");
                    let density = RadialTable::new(&r, &rho, &vec![1.0; r.len()]);
                    let (r_min, r_max) = (r[0], r[r.len() - 1]);
                    Ok(ClusterGas::Radial(hse::hydrostatic_table(|x| density.density(x), &pot, gamma, r_min, r_max, 4000)))
                } else {
                    Ok(ClusterGas::Radial(RadialTable::from_temperature(&r, &rho, &t, *mean_molecular_weight, gamma)))
                }
            }
        }
    }

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
        // Integrated orbits run a little past t_end so the table covers it after rounding.
        let integrated: Vec<(f64, f64)> = self
            .galaxy
            .iter()
            .filter_map(|g| match g.orbit {
                OrbitConfig::Integrate { dt_myr, output_dt_myr, .. } => Some((dt_myr * MYR, output_dt_myr * MYR)),
                _ => None,
            })
            .collect();
        let dt = integrated.iter().map(|x| x.0).fold(f64::INFINITY, f64::min);
        let output_dt = integrated.iter().map(|x| x.1).fold(f64::INFINITY, f64::min);
        let t_integrate = if integrated.is_empty() { t_end } else { t_end + output_dt };

        let mut motions = Vec::with_capacity(self.galaxy.len());
        for (k, g) in self.galaxy.iter().enumerate() {
            motions.push(match &g.orbit {
                OrbitConfig::ConstantVelocity { position_kpc: x, velocity_km_s: v } => {
                    let at = |t: f64| [x[0] + v[0] * t, x[1] + v[1] * t, x[2] + v[2] * t];
                    let t1 = (1.01 * t_end).max(t_integrate);
                    Motion::Prescribed(OrbitTable::new(vec![0.0, t1], vec![at(0.0), at(t1)], vec![*v, *v]).map_err(|e| ConfigError(e.to_string()))?)
                }
                OrbitConfig::Table { file } => {
                    let f = std::fs::File::open(file).map_err(|e| ConfigError(format!("{}: {e}", file.display())))?;
                    let table = OrbitTable::read_csv(std::io::BufReader::new(f)).map_err(|e| ConfigError(format!("{}: {e}", file.display())))?;
                    if table.t_start() > 0.0 || table.t_end() < t_end {
                        return err(format!(
                            "{}: orbit covers t = {:.1}-{:.1} Myr but the run needs 0-{:.1} Myr",
                            self.galaxy_name(k),
                            table.t_start() / MYR,
                            table.t_end() / MYR,
                            self.run.t_end_myr
                        ));
                    }
                    Motion::Prescribed(table)
                }
                OrbitConfig::Integrate { position_kpc, velocity_km_s, .. } => Motion::Integrate { x0: *position_kpc, v0: *velocity_km_s },
            });
        }
        if integrated.is_empty() {
            return Ok(motions.into_iter().map(|m| match m {
                Motion::Prescribed(t) => t,
                Motion::Integrate { .. } => unreachable!(),
            }).collect());
        }

        // Prescribed tables must reach the end of the integration; stop it early if one
        // does not (it still covers t_end, checked above).
        let t1 = motions
            .iter()
            .filter_map(|m| match m {
                Motion::Prescribed(t) => Some(t.t_end()),
                Motion::Integrate { .. } => None,
            })
            .fold(t_integrate, f64::min);
        let potentials: Vec<Composite> = (0..self.galaxy.len()).map(|k| self.galaxy_potential(k)).collect();
        let bodies: Vec<Body> = motions
            .into_iter()
            .enumerate()
            .map(|(k, motion)| Body { potential: &potentials[k], orientation: self.galaxy_orientation(k), motion })
            .collect();
        let cluster = self.cluster_potential();
        let cluster_ref = cluster.as_ref().map(|c| c as &dyn Potential);
        Ok(orbit::integrate_bodies(cluster_ref, &bodies, self.orbits.mutual_gravity, 0.0, t1, dt, output_dt.min(t1)))
    }
}

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
            if line.replace(' ', "") != "r_kpc,density_g_cm3,temperature_k" {
                return err(format!("{}: expected header r_kpc,density_g_cm3,temperature_k", path.display()));
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


#[cfg(test)]
mod tests {
    use super::*;

    const MINIMAL: &str = r#"
        [run]
        output_dir = "out"
        t_end_myr = 100

        [box]
        side_kpc = 80

        [cluster.gas]
        model = "uniform"
        density_g_cm3 = 1e-27
        temperature_k = 5e7

        [[galaxy]]
        name = "parent"

        [[galaxy.potential]]
        type = "hernquist"
        mass_msun = 1e11
        a_kpc = 1.0

        [galaxy.disk]
        mass_msun = 5e9
        scale_length_kpc = 4
        r_max_kpc = 20
        sound_speed_km_s = 20
        n_particles = 5000

        [galaxy.orbit]
        mode = "constant_velocity"
        position_kpc = [0, 0, 0]
        velocity_km_s = [0, 0, 1000]
    "#;

    /// A second galaxy appended to MINIMAL.
    const COMPANION: &str = r#"
        [[galaxy]]
        name = "companion"

        [[galaxy.potential]]
        type = "hernquist"
        mass_msun = 2e10
        a_kpc = 0.5

        [galaxy.disk]
        mass_msun = 1e9
        scale_length_kpc = 2
        r_max_kpc = 8
        sound_speed_km_s = 15
        n_particles = 300

        [galaxy.orbit]
        mode = "constant_velocity"
        position_kpc = [20, 0, 0]
        velocity_km_s = [0, 200, 1000]
    "#;

    fn parse(text: &str) -> Result<Config, ConfigError> {
        Config::from_str_in(text, Path::new("."))
    }

    #[test]
    fn minimal_config_parses_with_defaults() {
        let c = parse(MINIMAL).unwrap();
        assert_eq!(c.run.dt_max_myr, 2.0);
        assert_eq!(c.snapshot_every(), 5);
        assert_eq!(c.hydro.gamma, units::GAMMA);
        assert_eq!(c.rk4_steps_per_dt_max().unwrap(), 4);
        assert!(c.cluster_potential().is_none());
        assert!(c.ram_pressure().unwrap().is_none(), "ram pressure is off by default");
        let t = c.orbit_table(0).unwrap();
        assert!(t.t_end() >= c.t_end());
        assert_eq!(c.galaxy_name(0), "parent");
        assert_eq!(c.disk(0).unwrap().1, 5000);
        match c.cluster_gas().unwrap() {
            ClusterGas::Uniform { rho, .. } => assert!((rho * UNIT_DENSITY_CGS / 1e-27 - 1.0).abs() < 1e-12),
            _ => panic!("expected uniform gas"),
        }
    }

    #[test]
    fn unknown_keys_and_bad_values_are_errors() {
        let typo = MINIMAL.replace("side_kpc = 80", "sid_kpc = 80");
        let e = parse(&typo).unwrap_err();
        assert!(e.0.contains("sid_kpc") || e.0.contains("side_kpc"), "{e}");
        let unknown_type = MINIMAL.replace("type = \"hernquist\"", "type = \"plummer\"");
        assert!(parse(&unknown_type).is_err());
        let bad = MINIMAL.replace("t_end_myr = 100", "t_end_myr = -1");
        assert!(parse(&bad).unwrap_err().0.contains("t_end_myr"));
        let misaligned = MINIMAL.replace("[run]", "[run]\nsnapshot_interval_myr = 3");
        assert!(parse(&misaligned).unwrap_err().0.contains("multiple"));
        let beta_without_potential = MINIMAL.replace(
            "model = \"uniform\"\n        density_g_cm3 = 1e-27\n        temperature_k = 5e7",
            "model = \"beta\"\n        rho0_g_cm3 = 1e-26\n        r_c_kpc = 100\n        beta = 0.6",
        );
        assert!(parse(&beta_without_potential).unwrap_err().0.contains("cluster.potential"));
        let no_count = MINIMAL.replace("n_particles = 5000", "");
        assert!(parse(&no_count).unwrap_err().0.contains("n_particles"));
    }

    /// SPH-only sections of bigwig-vibe configs are not part of this schema.
    #[test]
    fn full_sph_sections_are_rejected() {
        for extra in ["[self_gravity]\nenabled = true\n", "[cooling]\nenabled = true\n", "[glass]\nn_particles = 100\n", "[feedback]\nenabled = true\n"] {
            assert!(parse(&format!("{MINIMAL}\n{extra}")).is_err(), "{extra}");
        }
        assert!(parse(&MINIMAL.replace("side_kpc = 80", "side_kpc = 80\n        recenter_fraction = 0.25")).is_err());
        assert!(parse(&MINIMAL.replace("[box]", "[hydro]\n        n_ngb = 64\n\n        [box]")).is_err());
    }

    #[test]
    fn time_dependent_parameters_and_unit_conversion() {
        let text = MINIMAL.replace(
            "mass_msun = 1e11\n        a_kpc = 1.0",
            "mass_msun = { t_myr = [0, 100], value = [1e11, 5e10] }\n        a_kpc = 1.0",
        );
        let c = parse(&text).unwrap();
        let p = c.galaxy_potential(0);
        let x = [50.0, 0.0, 0.0];
        let a0 = p.acceleration(&x, 0.0)[0];
        let a1 = p.acceleration(&x, 100.0 * MYR)[0];
        assert!((a1 / a0 - 0.5).abs() < 1e-12);
        let expected = Hernquist { mass: 10.0.into(), a: 1.0.into() }.acceleration(&x, 0.0)[0];
        assert!((a0 / expected - 1.0).abs() < 1e-12);
    }

    #[test]
    fn ram_pressure_and_galarp_sections() {
        let text = format!("{MINIMAL}\n[ram_pressure]\nenabled = true\nsurface_density_msun_pc2 = 20\n[ram_pressure.shadow]\nenabled = true\n[galarp]\ndt_myr = 0.25\n");
        let c = parse(&text).unwrap();
        let rp = c.ram_pressure().unwrap().unwrap();
        assert!(rp.shadow.is_some());
        assert_eq!(c.rk4_steps_per_dt_max().unwrap(), 8);
        assert!(parse(&text.replace("dt_myr = 0.25", "dt_myr = 0.3")).unwrap().rk4_steps_per_dt_max().is_err());
        assert!(parse(&text.replace("surface_density_msun_pc2 = 20", "surface_density_msun_pc2 = -1")).is_err());
        let none = text.replace("model = \"uniform\"\n        density_g_cm3 = 1e-27\n        temperature_k = 5e7", "model = \"none\"");
        assert!(parse(&none).unwrap_err().0.contains("ram_pressure"));
    }

    #[test]
    fn disk_options_from_the_disk_section() {
        let c = parse(MINIMAL).unwrap();
        let o = c.disk_options(0);
        assert!(!o.self_gravity && o.noise_fraction == 0.0);
        let sg = MINIMAL.replace("n_particles = 5000", "n_particles = 5000\n        include_self_gravity = true\n        noise_fraction = 0.1");
        let o2 = parse(&sg).unwrap().disk_options(0);
        assert!(o2.self_gravity);
        assert_eq!(o2.noise_fraction, 0.1);
    }

    #[test]
    fn galaxy_orientation_from_a_normal_or_angles() {
        let with = |extra: &str| MINIMAL.replacen("[[galaxy.potential]]", &format!("{extra}\n[[galaxy.potential]]"), 1);
        let c = parse(MINIMAL).unwrap();
        assert!(c.galaxy_orientation(0).is_identity());
        let n = parse(&with("disk_normal = [0, 1, 1]")).unwrap().galaxy_orientation(0).normal();
        assert!((n[1] - 0.5f64.sqrt()).abs() < 1e-12 && (n[2] - 0.5f64.sqrt()).abs() < 1e-12);
        let n = parse(&with("inclination_deg = 90\nazimuth_deg = 90")).unwrap().galaxy_orientation(0).normal();
        assert!((n[1] - 1.0).abs() < 1e-12);
        assert!(parse(&with("disk_normal = [0, 0, 1]\ninclination_deg = 10")).is_err(), "both");
        assert!(parse(&with("disk_normal = [0, 0, 0]")).is_err(), "zero");
        assert!(parse(&with("azimuth_deg = 10")).is_err(), "azimuth alone");
        assert!(parse(&with("inclination_deg = 200")).is_err(), "out of range");
    }

    #[test]
    fn no_icm_model_parses() {
        let text = MINIMAL.replace("model = \"uniform\"\n        density_g_cm3 = 1e-27\n        temperature_k = 5e7", "model = \"none\"");
        let c = parse(&text).unwrap();
        assert!(matches!(c.cluster_gas().unwrap(), ClusterGas::Uniform { rho, .. } if rho == 0.0));
    }

    #[test]
    fn a_second_galaxy_has_its_own_cloud_count() {
        let c = parse(&format!("{MINIMAL}{COMPANION}")).unwrap();
        assert_eq!(c.galaxy.len(), 2);
        assert_eq!(c.galaxy_name(1), "companion");
        assert_eq!(c.disk(1).unwrap().1, 300, "clouds need not share a mass");
        let p0 = c.orbit_table(0).unwrap().at(0.0).unwrap().0;
        let p1 = c.orbit_table(1).unwrap().at(0.0).unwrap().0;
        assert_eq!((p0[0], p1[0]), (0.0, 20.0));
        // Without any galaxy the config is rejected.
        let none = MINIMAL.split("[[galaxy]]").next().unwrap();
        assert!(parse(none).is_err());
    }

    #[test]
    fn integrated_galaxies_pull_on_each_other() {
        // Two equal Hernquist galaxies starting at rest 40 kpc apart, with no cluster: they
        // fall together symmetrically. Without mutual gravity nothing would pull them.
        let pair = |x: f64| {
            format!(
                "[[galaxy]]\nname = \"g{x}\"\n[[galaxy.potential]]\ntype = \"hernquist\"\nmass_msun = 1e11\na_kpc = 1.0\n[galaxy.orbit]\nmode = \"integrate\"\nposition_kpc = [{x}, 0, 0]\nvelocity_km_s = [0, 0, 0]\n"
            )
        };
        let base = MINIMAL.split("[[galaxy]]").next().unwrap().replace("t_end_myr = 100", "t_end_myr = 200");
        let c = parse(&format!("{base}{}{}", pair(20.0), pair(-20.0))).unwrap();
        let t = c.orbit_tables().unwrap();
        let (a, b) = (t[0].at(c.t_end()).unwrap().0, t[1].at(c.t_end()).unwrap().0);
        assert!(a[0] < 19.0 && (a[0] + b[0]).abs() < 1e-9, "{a:?} {b:?}");
        let off = format!("{base}[orbits]\nmutual_gravity = false\n{}{}", pair(20.0), pair(-20.0));
        assert!(parse(&off).unwrap_err().0.contains("cluster.potential"));
    }

    #[test]
    fn integrated_orbit_needs_and_uses_the_cluster_potential() {
        let text = MINIMAL
            .replace("mode = \"constant_velocity\"", "mode = \"integrate\"")
            .replace("[cluster.gas]", "[[cluster.potential]]\n        type = \"nfw_virial\"\n        m_vir_msun = 1e15\n        r_vir_kpc = 2000\n        concentration = 5\n\n        [cluster.gas]")
            .replace("position_kpc = [0, 0, 0]", "position_kpc = [1000, 0, 0]")
            .replace("velocity_km_s = [0, 0, 1000]", "velocity_km_s = [0, 300, 0]");
        let c = parse(&text).unwrap();
        let t = c.orbit_table(0).unwrap();
        let (x, _) = t.at(c.t_end()).unwrap();
        assert!(x[0] < 1000.0, "galaxy should fall inward");
    }

    #[test]
    fn typo_inside_a_potential_term_is_reported_with_its_location() {
        let text = MINIMAL.replace("a_kpc = 1.0", "a_kcp = 1.0");
        let e = parse(&text).unwrap_err();
        assert!(e.0.contains("a_kcp"), "{e}");
    }
}
