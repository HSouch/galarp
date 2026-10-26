use galarp_io::config::{ Config, ClusterGasConfig };


pub struct Preview {
    pub t_myr: Vec<f64>,

    pub pos: Vec<[f64; 3]>,
    pub p_ram_dyn_cm2: Vec<f64>,

    pub phi_dw_deg: Vec<f64>, 
    pub speed_km_s: Vec<f64>,
    pub n_h_cm3: Vec<f64>,
    /// Core radius of a beta-model ICM, drawn on the orbit.
    pub r_c_kpc: Option<f64>,
    pub ram_pressure_on: bool,
}

/// Number of samples: enough for smooth braille curves on a wide terminal.
const SAMPLES: usize = 400;

impl Preview {
    pub fn new(cfg: &Config) -> Result<Self, String> {
        let rows = galarp_io::preview::preview(cfg, SAMPLES).map_err(|e| e.0)?;
        let normal = galarp_io::preview::disk_normal(cfg, 0);
        let rows: Vec<_> = rows.into_iter().filter(|r| r.galaxy == 0).collect();
        let r_c_kpc = match cfg.cluster.gas {
            ClusterGasConfig::Beta { r_c_kpc, .. } => Some(r_c_kpc),
            _ => None,
        };
        Ok(Self {
            t_myr: rows.iter().map(|r| r.t_myr).collect(),
            pos: rows.iter().map(|r| r.pos).collect(),
            p_ram_dyn_cm2: rows.iter().map(|r| r.p_ram_dyn_cm2).collect(),
            phi_dw_deg: rows.iter().map(|r| phi_dw_deg(&r.vel, &normal)).collect(),
            speed_km_s: rows.iter().map(|r| r.speed_km_s).collect(),
            n_h_cm3: rows.iter().map(|r| r.n_h_icm_cm3).collect(),
            r_c_kpc,
            ram_pressure_on: cfg.ram_pressure.enabled,
        })
    }
}

/// phi_dw = arccos(-w . n) with w = -v/|v| the wind direction, i.e. arccos(v . n / |v|).
pub fn phi_dw_deg(vel: &[f64; 3], normal: &[f64; 3]) -> f64 {
    let speed = (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]).sqrt();
    if speed == 0.0 {
        return f64::NAN;
    }
    let cos = (vel[0] * normal[0] + vel[1] * normal[1] + vel[2] * normal[2]) / speed;
    cos.clamp(-1.0, 1.0).acos().to_degrees()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn phi_dw_is_signed() {
        let n = [0.0, 0.0, 1.0];
        // Moving along +n: the wind blows onto the +n face.
        assert_eq!(phi_dw_deg(&[0.0, 0.0, 1000.0], &n), 0.0);
        assert!((phi_dw_deg(&[1000.0, 0.0, 0.0], &n) - 90.0).abs() < 1e-12);
        assert_eq!(phi_dw_deg(&[0.0, 0.0, -5.0], &n), 180.0);
        assert!(phi_dw_deg(&[0.0; 3], &n).is_nan());
    }
}
