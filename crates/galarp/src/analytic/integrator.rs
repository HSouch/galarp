//! Runge-Kutta 4 integration
//! 


use rayon::prelude::*;

use super::clouds::Clouds;
use super::ram_pressure::RamPressure;
use crate::external::ExternalField;


pub struct CloudForces<'a> {
    pub field: &'a ExternalField,
    
    pub ram_pressure: Option<&'a RamPressure>,

}


impl CloudForces<'_> {

    /// Total accelerations at time `t` of clouds at `pos` moving at `vel`, and the
    /// ram-pressure parts and shadow factors (empty without ram pressure).
    pub fn evaluate(&self, clouds: &Clouds, t: f64, pos: &[[f64; 3]], vel: &[[f64; 3]]) -> (Vec<[f64; 3]>, Vec<[f64; 3]>, Vec<f64>) {
        let centers = self.field.galaxy_centers(t);
        let mut acc: Vec<[f64; 3]> = pos
            .par_iter()
            .map(|x| self.field.acceleration_with_centers(x, t, &centers))
            .collect();
        
        let Some(rp) = self.ram_pressure else { return (acc, Vec::new(), Vec::new()) };
        
        let (a_ram, s) = rp.evaluate(self.field, t, pos, vel, &clouds.mass, &clouds.galaxy);
        
        for (a, r) in acc.iter_mut().zip(&a_ram) {
            for k in 0..3 {
                a[k] += r[k];
            }
        }
        (acc, a_ram, s)
    }


    /// Stores RP accelerations and shadowing at time `t` in  clouds
    /// (for output).
    pub fn record(&self, clouds: &mut Clouds, t: f64) {
        let (_, a_ram, s) = self.evaluate(clouds, t, &clouds.pos, &clouds.vel);
        let n = clouds.len();
        clouds.a_ram = if a_ram.is_empty() { vec![[0.0; 3]; n] } else { a_ram };
        clouds.shadow = if s.is_empty() { vec![1.0; n] } else { s };
    }
}


/// Advances every cloud from `t` to `t + dt` with classical RK4.
pub fn rk4_step(clouds: &mut Clouds, forces: &CloudForces, t: f64, dt: f64) {
    let x0 = clouds.pos.clone();
    let v0 = clouds.vel.clone();
    let along = |x: &[[f64; 3]], d: &[[f64; 3]], h: f64| -> Vec<[f64; 3]> {
        x.par_iter().zip(d).map(|(a, b)| [a[0] + h * b[0], a[1] + h * b[1], a[2] + h * b[2]]).collect()
    };
    let (a1, ..) = forces.evaluate(clouds, t, &x0, &v0);
    let (x2, v2) = (along(&x0, &v0, 0.5 * dt), along(&v0, &a1, 0.5 * dt));
    let (a2, ..) = forces.evaluate(clouds, t + 0.5 * dt, &x2, &v2);
    let (x3, v3) = (along(&x0, &v2, 0.5 * dt), along(&v0, &a2, 0.5 * dt));
    let (a3, ..) = forces.evaluate(clouds, t + 0.5 * dt, &x3, &v3);
    let (x4, v4) = (along(&x0, &v3, dt), along(&v0, &a3, dt));
    let (a4, ..) = forces.evaluate(clouds, t + dt, &x4, &v4);
    let h = dt / 6.0;
    clouds.pos.par_iter_mut().enumerate().for_each(|(i, x)| {
        for k in 0..3 {
            x[k] = x0[i][k] + h * (v0[i][k] + 2.0 * v2[i][k] + 2.0 * v3[i][k] + v4[i][k]);
        }
    });
    clouds.vel.par_iter_mut().enumerate().for_each(|(i, v)| {
        for k in 0..3 {
            v[k] = v0[i][k] + h * (a1[i][k] + 2.0 * a2[i][k] + 2.0 * a3[i][k] + a4[i][k]);
        }
    });
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::cluster::ClusterGas;
    use crate::orbit::OrbitTable;
    use crate::potentials::{Hernquist, PointMass};
    use crate::analytic::DensityAt;
    use crate::units::G;
    use std::f64::consts::PI;

    #[test]
    fn a_cloud_in_a_constant_wind_follows_the_quadratic_drag_solution() {
        // No gravity; the galaxy moves at V along +z through uniform gas. A cloud starting
        // with the galaxy's velocity is slowed (cluster frame) as v(t) = V / (1 + k V t),
        // k = rho / Sigma.
        let v_gal = 1000.0;
        let orbit = OrbitTable::new(vec![0.0, 1.0], vec![[0.0; 3], [0.0, 0.0, v_gal]], vec![[0.0, 0.0, v_gal]; 2]).unwrap();
        let field = ExternalField::new().with_galaxy(PointMass { mass: 0.0.into() }, Some(orbit));
        let (rho, sigma) = (1e-6, 1e-4);
        let rp = RamPressure { 
            icm: ClusterGas::Uniform { rho, u: 1.0 }, 
            surface_density: sigma.into(), 
            density_at: DensityAt::Galaxy, 
            shadow: None };

        let forces = CloudForces { field: &field, ram_pressure: Some(&rp) };
        let mut clouds = Clouds::default();
        clouds.push(0, 0, [0.0; 3], [0.0, 0.0, v_gal], 1e-5);

        let (dt, n) = (1e-4, 1000);
        for s in 0..n {
            rk4_step(&mut clouds, &forces, s as f64 * dt, dt);
        }
        let t = n as f64 * dt;
        let k = rho / sigma;
        let v = v_gal / (1.0 + k * v_gal * t);
        let z = (1.0 + k * v_gal * t).ln() / k;
        assert!((clouds.vel[0][2] / v - 1.0).abs() < 1e-9, "{} vs {v}", clouds.vel[0][2]);
        assert!((clouds.pos[0][2] / z - 1.0).abs() < 1e-9);
    }

    #[test]
    fn without_wind_a_circular_orbit_closes() {
        let galaxy = Hernquist { mass: 10.0.into(), a: 1.0.into() };
        let field = ExternalField::new().with_galaxy(galaxy, None);
        let forces = CloudForces { field: &field, ram_pressure: None };
        let r = 5.0;
        let v_c = (G * 10.0 * r).sqrt() / (r + 1.0);
        let period = 2.0 * PI * r / v_c;
        let mut clouds = Clouds::default();
        clouds.push(0, 0, [r, 0.0, 0.0], [0.0, v_c, 0.0], 1e-5);
        let n = 400;
        for s in 0..n {
            rk4_step(&mut clouds, &forces, s as f64 * period / n as f64, period / n as f64);
        }
        let x = clouds.pos[0];
        assert!(((x[0] - r).powi(2) + x[1].powi(2)).sqrt() < 1e-5 * r, "{x:?}");
    }
}
