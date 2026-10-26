

use std::io::{self, Write};

use bigwig::units::{self, MYR, UNIT_LENGTH_CM};

use crate::config::{Config, ConfigError};


/// One sample of one galaxy.
#[derive(Clone, Debug, PartialEq)]
pub struct PreviewRow {
    pub galaxy: usize,
    pub t_myr: f64,
    /// Cluster-frame position (kpc) and velocity (km/s).
    pub pos: [f64; 3],
    pub vel: [f64; 3],
    /// Distance from the cluster center (kpc) and speed (km/s).
    pub r_kpc: f64,
    pub speed_km_s: f64,
    /// ICM density (g/cm^3) and hydrogen number density (cm^-3, X = 0.76) at the galaxy.
    pub rho_icm_g_cm3: f64,
    pub n_h_icm_cm3: f64,
    /// ICM sound speed (km/s) and the galaxy's Mach number through it.
    pub c_icm_km_s: f64,
    pub mach: f64,
    /// Ram pressure rho |v|^2, in dyn/cm^2 and as P/k_B in K/cm^3.
    pub p_ram_dyn_cm2: f64,
    pub p_ram_over_k: f64,
    /// Angle between the wind and the disk normal, in [0, 90] degrees (90 = edge-on); NaN if
    /// the galaxy is at rest.
    pub wind_angle_deg: f64,
}

/// The disk's unit normal (spin axis) in the cluster frame.
pub fn disk_normal(cfg: &Config, k: usize) -> [f64; 3] {
    cfg.galaxy_orientation(k).normal()
}

/// Samples every galaxy at `n` evenly spaced times over [0, t_end].
pub fn preview(cfg: &Config, n: usize) -> Result<Vec<PreviewRow>, ConfigError> {
    let gas = cfg.cluster_gas()?;
    let t_end = cfg.t_end();
    let n = n.max(2);
    let mut rows = Vec::with_capacity(n * cfg.galaxy.len());
    for (k, orbit) in cfg.orbit_tables()?.into_iter().enumerate() {
        let normal = disk_normal(cfg, k);
        for s in 0..n {
            let t = t_end * s as f64 / (n - 1) as f64;
            let (pos, vel) = orbit.at(t).expect("orbit_table covers [0, t_end]");
            let r = norm(&pos);
            let speed = norm(&vel);
            let rho = gas.density_at(&pos);
            let u = gas.specific_energy_at(&pos);
            let c = (units::GAMMA * (units::GAMMA - 1.0) * u).sqrt();
            let rho_cgs = rho * units::UNIT_DENSITY_CGS;
            let p_ram = rho_cgs * (speed * units::UNIT_VELOCITY_CM_S).powi(2);
            // No motion, no wind: the angle is undefined (NaN).
            let cos = if speed > 0.0 { (dot(&vel, &normal) / speed).abs().min(1.0) } else { f64::NAN };
            rows.push(PreviewRow {
                galaxy: k,
                t_myr: t / MYR,
                pos,
                vel,
                r_kpc: r,
                speed_km_s: speed,
                rho_icm_g_cm3: rho_cgs,
                n_h_icm_cm3: 0.76 * rho_cgs / units::PROTON_MASS_G,
                c_icm_km_s: c,
                mach: if c > 0.0 { speed / c } else { f64::INFINITY },
                p_ram_dyn_cm2: p_ram,
                p_ram_over_k: p_ram / units::BOLTZMANN_CGS,
                wind_angle_deg: cos.acos().to_degrees(),
            });
        }
    }
    Ok(rows)
}

/// Writes the rows as CSV, with the galaxy names and box size in `#` comment lines.
pub fn write_preview_csv(cfg: &Config, rows: &[PreviewRow], mut w: impl Write) -> io::Result<()> {
    for k in 0..cfg.galaxy.len() {
        writeln!(w, "# galaxy {k}: {}", cfg.galaxy_name(k))?;
    }
    writeln!(w, "# box_side_kpc: {}", cfg.box_.side_kpc)?;
    writeln!(w, "# unit_length_cm: {UNIT_LENGTH_CM:e}")?;
    writeln!(
        w,
        "galaxy,t_myr,x_kpc,y_kpc,z_kpc,vx_km_s,vy_km_s,vz_km_s,r_kpc,speed_km_s,rho_icm_g_cm3,n_h_icm_cm3,c_icm_km_s,mach,p_ram_dyn_cm2,p_ram_over_k,wind_angle_deg"
    )?;
    for r in rows {
        writeln!(
            w,
            "{},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e},{:.6e}",
            r.galaxy,
            r.t_myr,
            r.pos[0],
            r.pos[1],
            r.pos[2],
            r.vel[0],
            r.vel[1],
            r.vel[2],
            r.r_kpc,
            r.speed_km_s,
            r.rho_icm_g_cm3,
            r.n_h_icm_cm3,
            r.c_icm_km_s,
            r.mach,
            r.p_ram_dyn_cm2,
            r.p_ram_over_k,
            r.wind_angle_deg
        )?;
    }
    Ok(())
}

fn norm(a: &[f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[cfg(test)]
mod tests {
    use super::*;

    const FACE_ON: &str = r#"
[run]
output_dir = "x"
t_end_myr = 100
dt_max_myr = 2
snapshot_interval_myr = 10
[box]
side_kpc = 80
[cluster.gas]
model = "uniform"
density_g_cm3 = 1e-27
temperature_k = 5e7
[[galaxy]]
name = "g"
[[galaxy.potential]]
type = "hernquist"
mass_msun = 1e10
a_kpc = 0.5
[galaxy.orbit]
mode = "constant_velocity"
position_kpc = [0, 0, 0]
velocity_km_s = [0, 600, 800]
"#;

    #[test]
    fn uniform_wind_gives_constant_ram_pressure_angle_and_mach() {
        let cfg = Config::from_str_in(FACE_ON, std::path::Path::new(".")).unwrap();
        let rows = preview(&cfg, 11).unwrap();
        assert_eq!(rows.len(), 11);
        let last = rows.last().unwrap();
        assert!((last.t_myr - 100.0).abs() < 1e-9);
        // 1000 km/s through 1e-27 g/cm^3: 1e-11 dyn/cm^2.
        assert!(rows.iter().all(|r| (r.p_ram_dyn_cm2 / 1e-11 - 1.0).abs() < 1e-9 && (r.speed_km_s - 1000.0).abs() < 1e-9));
        // v along (0, 0.6, 0.8): cos(angle to z) = 0.8, angle 36.87 degrees.
        assert!(rows.iter().all(|r| (r.wind_angle_deg - 0.8f64.acos().to_degrees()).abs() < 1e-9));
        // 5e7 K at mu = 0.6: c = sqrt(gamma k T / (mu m_p)) = 1072 km/s.
        assert!((last.c_icm_km_s / 1072.0 - 1.0).abs() < 2e-3, "{}", last.c_icm_km_s);
        // The galaxy moves 1000 km/s * 100 Myr = 102.3 kpc.
        assert!((last.r_kpc - 102.27).abs() < 0.01, "{}", last.r_kpc);
        let mut csv = Vec::new();
        write_preview_csv(&cfg, &rows, &mut csv).unwrap();
        let text = String::from_utf8(csv).unwrap();
        assert!(text.starts_with("# galaxy 0: g\n") && text.lines().count() == 3 + 1 + 11);
    }
}
