//! Radial profile table of density and specific internal energy.
//! interpolated linearly in log r, log rho and log u. 
//! Outside the table the end segments are extended as power laws.


use crate::units;


pub struct RadialTable {
    log_r: Vec<f64>,
    log_rho: Vec<f64>,
    log_u: Vec<f64>,
}

impl RadialTable {

    pub fn new(r: &[f64], rho: &[f64], u: &[f64]) -> Self {
        assert!(r.len() >= 2 && rho.len() == r.len() && u.len() == r.len(), "need matching columns of length >= 2");
        assert!(r.windows(2).all(|w| w[0] < w[1]), "radii must increase");
        assert!(r.iter().chain(rho).chain(u).all(|&v| v > 0.0), "values must be positive");
        Self {
            log_r: r.iter().map(|v| v.ln()).collect(),
            log_rho: rho.iter().map(|v| v.ln()).collect(),
            log_u: u.iter().map(|v| v.ln()).collect(),
        }
    }

    /// From temperatures in kelvin, for observed profiles.
    pub fn from_temperature(r: &[f64], rho: &[f64], t_kelvin: &[f64], mu: f64, gamma: f64) -> Self {
        let u: Vec<f64> = t_kelvin.iter().map(|&t| units::u_from_temperature(t, mu, gamma)).collect();
        Self::new(r, rho, &u)
    }

    pub fn r_min(&self) -> f64 {
        self.log_r[0].exp()
    }

    pub fn r_max(&self) -> f64 {
        self.log_r[self.log_r.len() - 1].exp()
    }

    pub fn density(&self, r: f64) -> f64 {
        self.interp(&self.log_rho, r)
    }

    pub fn specific_energy(&self, r: f64) -> f64 {
        self.interp(&self.log_u, r)
    }

    fn interp(&self, col: &[f64], r: f64) -> f64 {
        let lr = r.max(f64::MIN_POSITIVE).ln();
        let n = self.log_r.len();
        let k = self.log_r.partition_point(|&x| x <= lr).clamp(1, n - 1);
        let f = (lr - self.log_r[k - 1]) / (self.log_r[k] - self.log_r[k - 1]);
        (col[k - 1] + f * (col[k] - col[k - 1])).exp()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn power_laws_are_exact_inside_and_outside_the_table() {
        let r = [1.0, 2.0, 4.0, 8.0];
        let rho: Vec<f64> = r.iter().map(|x: &f64| 3.0 * x.powf(-2.0)).collect();
        let u: Vec<f64> = r.iter().map(|x: &f64| 5.0 * x.powf(0.5)).collect();
        let t = RadialTable::new(&r, &rho, &u);
        for x in [0.5, 1.0, 3.0, 8.0, 20.0] {
            assert!((t.density(x) / (3.0 * x.powf(-2.0)) - 1.0).abs() < 1e-12);
            assert!((t.specific_energy(x) / (5.0 * x.powf(0.5)) - 1.0).abs() < 1e-12);
        }
    }
}
