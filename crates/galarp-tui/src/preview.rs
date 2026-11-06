
use galarp_io::config::{Config, PotentialTerm};


/// Number of samples for generating curves 
const NSAMPLES: usize = 400;

pub struct Preview {
    pub t_myr: Vec<f64>,

    pub pos: Vec<[f64; 3]>,

    pub p_ram_dyn_cm2: Vec<f64>,

    pub phi_dw_deg: Vec<f64>,
    pub speed_km_s: Vec<f64>,
    pub n_h_cm3: Vec<f64>,

    pub r_vir_kpc: Option<f64>,
    pub ram_pressure_on: bool,

}


impl Preview {
    pub fn new(cfg: &Config) -> Result<Self, String> {
        let rows = galarp_io::preview::preview(cfg, NSAMPLES).map_err(|e| e.0)?;
        let normal = galarp_io::preview::disk_normal(cfg, 0);

        let rows: Vec<_> = rows.into_iter().filter(|r| r.galaxy == 0).collect();

        let r_vir_kpc = cfg.cluster.potential.iter()
            .find_map(|term| match term {
                PotentialTerm::NfwVirial {r_vir_kpc, .. } => Some(*r_vir_kpc),
                _ => None,
            });
        

        Ok(Self {
            t_myr: rows.iter().map(|row| row.t_myr).collect(),
            pos: rows.iter().map(|r| r.pos).collect(),
            p_ram_dyn_cm2: rows.iter().map(|row| row.p_ram_dyn_cm2).collect(),
            phi_dw_deg: rows.iter().map(|row| phi_dw_deg(&row.vel, &normal)).collect(),
            speed_km_s: rows.iter().map(|row| row.speed_km_s).collect(),
            n_h_cm3: rows.iter().map(|row| row.n_h_icm_cm3).collect(),
            r_vir_kpc,
            ram_pressure_on: cfg.ram_pressure.enabled,
        })
    }
}


pub fn phi_dw_deg(vel: &[f64; 3], normal: &[f64; 3]) -> f64 {
    let speed = (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]).sqrt();
    if speed == 0.0 {
        return f64::NAN;
    }

    let cos = (vel[0] * normal[0] + vel[1] * normal[1] + vel[2] * normal[2]) / speed;
    cos.clamp(-1.0, 1.0).acos().to_degrees()
}