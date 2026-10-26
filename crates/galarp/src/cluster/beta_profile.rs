//! Beta profile for ICM density
//! 


#[derive(Clone, Copy, Debug)]
pub struct BetaModel {
    pub rho_0: f64,         // Central density
    pub r_c: f64,           // Intrinsic size
    pub beta:  f64,         // Slope param
}

impl BetaModel {
    pub fn density(&self, r: f64) -> f64 {
        let s = r / self.r_c;
        self.rho_0 * (1.0 + s * s).powf(-1.5 * self.beta)
    }
}


#[cfg(test)]
pub mod tests{
    use super::*;

    #[test]
    fn profile_returns_appropriate_central_density() {
        let rho_0 = 1e-5;
        let model = BetaModel{rho_0: rho_0, r_c: 1.0, beta: 0.5};
        assert!((model.density(0.0) - rho_0) / rho_0 < 1e-6);
    }

    #[test]
    fn profile_close_to_zero_at_large_r() {
        let model = BetaModel{rho_0: 1.0, r_c: 1.0, beta: 0.5};
        let dens = model.density(1e5);
        assert!(dens < 1e-6, "dens = {dens}");
    }
}