
// CODE UNITS

pub const UNIT_LENGTH_CM: f64 = 3.085_677_581e21;   // Code length unit: 1kpc
pub const UNIT_MASS_G: f64 = 1.988_92e43;           // Code mass unit: 1 Msun
pub const UNIT_VELOCITY_CM_S: f64 = 1.0e5;          // Code velocity unit 1 km/s
pub const UNIT_TIME_S: f64 = UNIT_LENGTH_CM / UNIT_VELOCITY_CM_S;

pub const UNIT_SURFACE_DENSITY_CGS: f64 = UNIT_MASS_G / (UNIT_LENGTH_CM * UNIT_LENGTH_CM);
pub const UNIT_DENSITY_CGS: f64 = UNIT_MASS_G / (UNIT_LENGTH_CM * UNIT_LENGTH_CM * UNIT_LENGTH_CM);

pub const UNIT_SPECIFIC_ENERGY_CGS: f64 = UNIT_VELOCITY_CM_S * UNIT_VELOCITY_CM_S;


// PHYSICAL CONSTANTS

pub const G_CGS: f64 = 6.674_30e-8;
pub const BOLTZMANN_CGS: f64 = 1.380_649e-16;
pub const PROTON_MASS_G: f64 = 1.672_621_924e-24;

pub const G: f64 = G_CGS * UNIT_MASS_G / (UNIT_LENGTH_CM * UNIT_VELOCITY_CM_S * UNIT_VELOCITY_CM_S);

pub const SOLAR_METALLICITY: f64 = 0.0134;

pub const GAMMA: f64 = 5.0 / 3.0;                   // Adiabatic constant for a monatomic ideal gas


/// Determines specific energy from a given temperature with mean molecular weight mu
pub fn u_from_temperature(t_kelvlin: f64, mu: f64, gamma: f64) -> f64 {
    BOLTZMANN_CGS * t_kelvlin / ((gamma - 1.0) * mu * PROTON_MASS_G) / UNIT_SPECIFIC_ENERGY_CGS
}



/// Determines temperature in Kelvin of gas with specific internal energy u (code units).
pub fn temperature_from_u(u: f64, mu: f64, gamma: f64) -> f64 {
    u * UNIT_SPECIFIC_ENERGY_CGS * (gamma - 1.0) * mu * PROTON_MASS_G / BOLTZMANN_CGS
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_against_astropy() {
        assert!((G_CGS - 6.6743e-8).abs() / 6.6743e-8 < 1e-3, "G_CGS = {G_CGS}");
        assert!((BOLTZMANN_CGS - 1.380649e-16).abs() / 1.380649e-16 < 1e-3, "k_B = {BOLTZMANN_CGS}");
        assert!((PROTON_MASS_G - 1.6726e-24).abs() / 1.6726e-24 < 1e-3, "mu_p = {PROTON_MASS_G}");
    }

    #[test]
    fn grav_constant_matches_gadget() {
        assert!((G - 43_007.1).abs() / 43_007.1 < 1e-3, "G = {G}")
    }

    #[test]
    fn temp_round_trip() {
        let u = u_from_temperature(1e7, 0.6, GAMMA);
        assert!((temperature_from_u(u, 0.6, GAMMA) - 1e7).abs() < 1e-6);
    }
}