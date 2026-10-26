use std::ops::ControlFlow;
use std::path::Path;

use bigwig::external::{ExternalField, GalaxyField};
use bigwig::ic::GasDisk;
use bigwig::galarp::{rk4_step, CloudForces, Clouds, RamPressure};

use bigwig::units::{self, MYR};
use bigwig_io::config::{Config};
use bigwig_io::log::CsvLog;
use bigwig_io::snapshot::{output_path, write_cloud_snapshot, SnapshotHeader, WriteOptions};


pub type Error = Box<dyn std::error::Error + Send + Sync>;

pub const GALAXY_COLUMNS: [&str; 7] = [
    "in_disk_fraction",
    "bound_fraction",
    "mean_shadow",
    "rho_icm_g_cm3",
    "wind_km_s",
    "p_ram_dyn_cm2",
    "wind_angle_deg",
];


pub struct Run {
    pub cfg: Config,
    config_text: String,
    field: ExternalField,
    rp: Option<RamPressure>,
    pub clouds: Clouds,

    pub summary: Vec<String>,
    steps_per_chunk: u64,
}



impl Run {
    pub fn new(cfg: Config, config_text: String) -> Result<Self, Error> {
        let steps_per_chunk = cfg.rk4_steps_per_dt_max()?;
        let orbits = cfg.orbit_tables()?;

        let mut field = ExternalField::new();

        for (k, orbit) in orbits.into_iter().enumerate() {
            field.galaxies.push(
                GalaxyField::new(cfg.galaxy_potential(k),
                Some(orbit)).with_orientation(cfg.galaxy_orientation(k)));
        }
        
        if let Some(p) = cfg.cluster_potential() {
            field.cluster = Some(Box::new(p));
        }
        let rp = cfg.ram_pressure()?;

        let mut summary = Vec::new();

        let mut clouds = Clouds::default();

        let mut next_id = 0;

        for k in 0..cfg.galaxy.len() {
            let Some((params, n)) = cfg.disk(k) else {
                summary.push(format!("  {}: no disk (gravity only)", cfg.galaxy_name(k)));
                continue;
            };
            let potential = cfg.galaxy_potential(k);
            let disk = GasDisk::with_options(params, &potential, cfg.disk_options(k));
            let (x, v) = field.galaxies[k].state(0.0);
            let sigma_z = cfg.galarp.vertical_dispersion_km_s.unwrap_or(params.sound_speed);
            let seed = cfg.galarp.seed.wrapping_add(k as u64);
            clouds.add_disk(
                &disk, 
                &potential, 
                cfg.galaxy_orientation(k), 
                k as u16, 
                n, 
                x, 
                v, 
                sigma_z, 
                seed, 
                &mut next_id);
            summary.push(
                format!("  {}: {n} clouds of {:.3e} Msun, sigma_z {sigma_z} km/s",
                cfg.galaxy_name(k), params.mass / n as f64 / units::MSUN));
        }   
        if clouds.is_empty() {
            return Err("no galaxy has a [galaxy.disk], so no clouds".into());
        }
        
        let  rpc = &cfg.ram_pressure;

        summary.push(match &rp {
            Some(_) => format!(
                "ram pressure: Sigma {} Msun/pc^2, rho at the {}, shadow {}",
                rpc.surface_density_msun_pc2,
                rpc.density_at,
                if rpc.shadow.enabled { format!("C_D {} in {} kpc columns", rpc.shadow.drag_coefficient, rpc.shadow.bin_width_kpc) } else { "off".into() }
            ),
            None => "ram pressure: off (gravity only)".into(),
        });

        summary.push(format!(
            "run {} Myr in RK4 steps of {} Myr, log every {} Myr, snapshot every {} Myr",
            cfg.run.t_end_myr, cfg.galarp.dt_myr, cfg.run.dt_max_myr, cfg.run.snapshot_interval_myr
        ));
        Ok(Self { cfg, config_text, field, rp, clouds, summary, steps_per_chunk })
    }

    /// Number of log steps after the initial state.
    pub fn chunks(&self) -> u64 {
        (self.cfg.t_end() / self.cfg.dt_max()).round() as u64
    }

    pub fn log_header(&self) -> Vec<String> {
        let mut h = vec!["t_myr".to_string()];
        for k in 0..self.cfg.galaxy.len() {
            h.extend(GALAXY_COLUMNS.map(|c| format!("g{k}_{c}")));
        }
        h
    }

    /// Integrates to `run.t_end_myr`, writing into `run.output_dir` (which must be free of
    /// snapshots or be meant to be overwritten), and calls `observe` after the initial state
    /// and every log step. Returns false if `observe` stopped the run early.
   
    pub fn execute(&mut self, mut observe: impl FnMut(&Progress) -> ControlFlow<()>) -> Result<bool, Error> {
        // Handle outdir and write config toml file
        let out = self.cfg.run.output_dir.clone();
        std::fs::create_dir_all(&out)?;
        std::fs::write(out.join("config.toml"), &self.config_text)?;
        
        // Create CSV log
        let header = self.log_header();
        let mut log = CsvLog::create(out.join("log.csv"), &header.iter().map(String::as_str).collect::<Vec<_>>())?;
        
        
        let Self { cfg, config_text, field, rp, clouds, steps_per_chunk, .. } = self;

        let forces = CloudForces { field, ram_pressure: rp.as_ref() };
        let dt = cfg.galarp.dt_myr * MYR;
        let chunks = (cfg.t_end() / cfg.dt_max()).round() as u64;
        let snap_every = cfg.snapshot_every();
        for chunk in 0..=chunks {
            if chunk > 0 {
                for s in 0..*steps_per_chunk {
                    let step = (chunk - 1) * *steps_per_chunk + s;
                    rk4_step(clouds, &forces, step as f64 * dt, dt);
                }
            }
            let t = (chunk * *steps_per_chunk) as f64 * dt;
            forces.record(clouds, t);
            let row = log_row(cfg, field, rp.as_ref(), clouds, t);
            log.row(&row)?;
            let states: Vec<_> = field.galaxies.iter().map(|g| g.state(t)).collect();
            let galaxy_coords: Vec<[f64; 3]> = (0..clouds.len()).map(|i| galaxy_frame(field, &states, clouds.galaxy[i] as usize, &clouds.pos[i])).collect();
            let snapshot = (chunk % snap_every == 0).then_some(chunk / snap_every);
            if let Some(number) = snapshot {
                write_output(cfg, config_text, &out, field, &states, clouds, &galaxy_coords, number, t)?;
            }
            let wind_dir: Vec<Option<[f64; 3]>> = field
                .galaxies
                .iter()
                .zip(&states)
                .map(|(g, (_, v))| {
                    let s = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
                    (rp.is_some() && s > 0.0).then(|| g.orientation.to_galaxy(&[-v[0] / s, -v[1] / s, -v[2] / s]))
                })
                .collect();
            let p = Progress { 
                chunk, 
                chunks, 
                t_myr: t / MYR, 
                row: &row, 
                snapshot, 
                clouds, 
                galaxy_coords: &galaxy_coords, 
                wind_dir: &wind_dir 
            };
            if observe(&p).is_break() {
                return Ok(false);
            }
        }
        Ok(true)
    }
}



pub struct Progress<'a> {
    pub chunk: u64,
    pub chunks: u64,

    pub t_myr: f64,
    pub row: &'a [f64],

    pub snapshot: Option<u64>,
    pub clouds: &'a Clouds,

    pub galaxy_coords: &'a [[f64; 3]],

    pub wind_dir: &'a [Option<[f64; 3]>],
}


#[allow(clippy::too_many_arguments)]
fn write_output(
    cfg: &Config,
    config_text: &str,
    out: &Path,
    field: &ExternalField,
    states: &[([f64; 3], [f64; 3])],
    clouds: &Clouds,
    galaxy_coords: &[[f64; 3]],
    number: u64,
    t: f64,
) -> Result<(), Error> {
    let side = cfg.box_.side_kpc;
    let parent = states[0].0;
    let header = SnapshotHeader {
        time: t,
        box_size: side,
        box_origin: [parent[0] - 0.5 * side, parent[1] - 0.5 * side, parent[2] - 0.5 * side],
        galaxy_pos: states.iter().map(|x| x.0).collect(),
        galaxy_vel: states.iter().map(|x| x.1).collect(),
        galaxy_normal: field.galaxies.iter().map(|g| g.orientation.normal()).collect(),
        gamma: cfg.hydro.gamma,
        mean_molecular_weight: cfg.run.mean_molecular_weight,
    };
    let path = output_path(out, number);
    std::fs::create_dir_all(path.parent().unwrap())?;
    write_cloud_snapshot(&path, &header, clouds, galaxy_coords, config_text, WriteOptions { compress: cfg.run.compress })?;
    Ok(())
}

/// Offset of cluster position `x` from galaxy `g`'s center, in that galaxy's frame.
fn galaxy_frame(field: &ExternalField, states: &[([f64; 3], [f64; 3])], g: usize, x: &[f64; 3]) -> [f64; 3] {
    let c = states[g].0;
    field.galaxies[g].orientation.to_galaxy(&[x[0] - c[0], x[1] - c[1], x[2] - c[2]])
}

/// Mass fraction of each galaxy's clouds bound to it: 0.5 |v - v_galaxy|^2 plus the
/// galaxy's own potential (not the cluster's) below zero.
fn bound_fraction(field: &ExternalField, clouds: &Clouds, t: f64) -> Vec<f64> {
    let states: Vec<_> = field.galaxies.iter().map(|g| g.state(t)).collect();
    let n_gal = field.galaxies.len();
    let (mut bound, mut total) = (vec![0.0; n_gal], vec![0.0; n_gal]);
    for i in 0..clouds.len() {
        let g = clouds.galaxy[i] as usize;
        let (c, vg) = states[g];
        let (x, v) = (clouds.pos[i], clouds.vel[i]);
        let dv = [v[0] - vg[0], v[1] - vg[1], v[2] - vg[2]];
        let e = 0.5 * (dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]) + field.galaxies[g].potential_at_offset(&[x[0] - c[0], x[1] - c[1], x[2] - c[2]], t);
        total[g] += clouds.mass[i];
        if e < 0.0 {
            bound[g] += clouds.mass[i];
        }
    }
    bound.iter().zip(&total).map(|(a, b)| if *b > 0.0 { a / b } else { f64::NAN }).collect()
}

/// Mass fraction of each galaxy's clouds inside its `[diagnostics]` disk region.
fn in_disk_fraction(cfg: &Config, field: &ExternalField, clouds: &Clouds, t: f64) -> Vec<f64> {
    let (r_reg, z_reg) = (cfg.diagnostics.disk_region_radius_kpc, cfg.diagnostics.disk_region_half_height_kpc);
    let states: Vec<_> = field.galaxies.iter().map(|g| g.state(t)).collect();
    let n_gal = cfg.galaxy.len();
    let (mut inside, mut total) = (vec![0.0; n_gal], vec![0.0; n_gal]);
    for i in 0..clouds.len() {
        let g = clouds.galaxy[i] as usize;
        let d = galaxy_frame(field, &states, g, &clouds.pos[i]);
        total[g] += clouds.mass[i];
        if d[0].hypot(d[1]) < r_reg && d[2].abs() < z_reg {
            inside[g] += clouds.mass[i];
        }
    }
    inside.iter().zip(&total).map(|(a, b)| if *b > 0.0 { a / b } else { f64::NAN }).collect()
}

fn log_row(cfg: &Config, field: &ExternalField, rp: Option<&RamPressure>, clouds: &Clouds, t: f64) -> Vec<f64> {
    let frac = in_disk_fraction(cfg, field, clouds, t);
    let bound = bound_fraction(field, clouds, t);
    let mut row = vec![t / MYR];
    for (k, g) in field.galaxies.iter().enumerate() {
        let (mut sm, mut m) = (0.0, 0.0);
        for i in (0..clouds.len()).filter(|&i| clouds.galaxy[i] as usize == k) {
            sm += clouds.shadow[i] * clouds.mass[i];
            m += clouds.mass[i];
        }
        let (center, v) = g.state(t);
        let speed = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        let rho = rp.map_or(0.0, |rp| rp.icm.density_at(&center));
        let angle = match rp {
            Some(rp) => rp.winds(field, t)[k].angle_deg(),
            None => f64::NAN,
        };
        row.extend([
            frac[k], 
            bound[k], if m > 0.0 { sm / m } else { f64::NAN }, 
            rho * units::UNIT_DENSITY_CGS, speed, rho * speed * speed * units::UNIT_PRESSURE_CGS, 
            angle
            ]);
    }
    row
}
