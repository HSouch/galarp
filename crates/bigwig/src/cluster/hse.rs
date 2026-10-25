
//! Temperature from hydrostatic equilibrium in the cluster potential.
//!


use super::RadialTable;
use crate::potentials::Potential;

/// Inward gravitational acceleration at radius `r`, measured along +x at t = 0.
pub fn inward_gravity(potential: &dyn Potential, r: f64) -> f64 {
    -potential.acceleration(&[r, 0.0, 0.0], 0.0)[0]
}

/// Builds a radial table of `density(r)` and the 
/// hydrostatic specific energy on `n`
/// log-spaced radii from `r_min` to `r_max`.
pub fn hydrostatic_table(
    density: impl Fn(f64) -> f64,
    potential: &dyn Potential,
    gamma: f64,
    r_min: f64,
    r_max: f64,
    n: usize,
) -> RadialTable {
    assert!(n >= 2 && r_min > 0.0 && r_max > r_min);
    let dlog = (r_max / r_min).ln() / (n - 1) as f64;
    let r: Vec<f64> = (0..n).map(|i| r_min * (i as f64 * dlog).exp()).collect();
    let rho: Vec<f64> = r.iter().map(|&x| density(x)).collect();
    // Integrand in ln r: rho g r.
    let f: Vec<f64> = (0..n).map(|i| rho[i] * inward_gravity(potential, r[i]) * r[i]).collect();
    let mut p = vec![0.0; n];
    for i in (0..n - 1).rev() {
        p[i] = p[i + 1] + 0.5 * (f[i] + f[i + 1]) * dlog;
    }
    // The outermost point has P = 0 by construction; give it the next point's u so the
    // table stays positive (it only matters beyond r_max, which callers should avoid).
    let mut u: Vec<f64> = (0..n).map(|i| p[i] / ((gamma - 1.0) * rho[i])).collect();
    u[n - 1] = u[n - 2];
    RadialTable::new(&r, &rho, &u)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cluster::{BetaModel, ClusterGas};
    use crate::potentials::{Nfw, PointMass};
    use crate::units::{G, GAMMA};

    #[test]
    fn singular_isothermal_density_around_a_point_mass_matches_the_exact_solution() {
        // rho = A / r^2 and g = G M / r^2 give P = A G M / (3 r^3), so u = G M / (3 (gamma - 1) r).
        let p = PointMass { mass: 2.0.into() };
        let table = hydrostatic_table(|r| 5.0 / (r * r), &p, GAMMA, 1.0, 1e4, 4000);
        for r in [1.0, 3.0, 10.0, 100.0] {
            let exact = G * 2.0 / (3.0 * (GAMMA - 1.0) * r);
            assert!((table.specific_energy(r) / exact - 1.0).abs() < 1e-4, "r = {r}");
        }
    }

    #[test]
    fn beta_model_in_nfw_cluster_satisfies_hydrostatic_equilibrium() {
        let cluster = Nfw::from_virial(1e5, 2000.0, 5.0);
        let beta = BetaModel { rho_0: 1e-3, r_c: 100.0, beta: 0.6 };
        let gas = ClusterGas::beta_model_hse(&beta, &cluster, GAMMA, 1.0, 1e5, 4000);
        let pressure = |r: f64| (GAMMA - 1.0) * gas.density(r) * gas.specific_energy(r);
        for r in [10.0, 100.0, 500.0, 2000.0] {
            let dr = 1e-3 * r;
            let dp_dr = (pressure(r + dr) - pressure(r - dr)) / (2.0 * dr);
            let expected = -gas.density(r) * inward_gravity(&cluster, r);
            assert!((dp_dr / expected - 1.0).abs() < 1e-3, "r = {r}: {dp_dr} vs {expected}");
            assert!((gas.density(r) / beta.density(r) - 1.0).abs() < 1e-6);
        }
        // A 1e15 Msun cluster should be several keV inside the core.
        let t = crate::units::temperature_from_u(gas.specific_energy(100.0), 0.6, GAMMA);
        assert!(t > 3e7 && t < 3e8, "T(100 kpc) = {t:.3e} K");
    }
}
