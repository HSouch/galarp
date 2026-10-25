//! Integrator for galaRP
//! 
//! We are using a simple leapfrog integrator for energy preservation over many orbits
//! 


use super::table::OrbitTable;
use crate::potentials::Potential;


pub fn integrate(
    potential: &dyn Potential,
    x0: [f64; 3],
    v0: [f64; 3],
    t0: f64,
    t_end: f64,
    dt: f64,
    output_dt: f64,
) -> OrbitTable {
    assert!(t_end > t0 && dt > 0.0 && output_dt > 0.0, "need t_end > t0 and positive steps");
    let n_out = ((t_end - t0) / output_dt).round().max(1.0) as usize;
    let out_dt = (t_end - t0) / n_out as f64;
    let n_sub = (out_dt / dt).ceil().max(1.0) as usize;
    let h = out_dt / n_sub as f64;

    let (mut x, mut v) = (x0, v0);
    let mut t = t0;
    let mut a = potential.acceleration(&x, t);
    let (mut ts, mut xs, mut vs) = (vec![t], vec![x], vec![v]);
    for k in 1..=n_out {
        for _ in 0..n_sub {
            for d in 0..3 {
                v[d] += 0.5 * h * a[d];
                x[d] += h * v[d];
            }
            t += h;
            a = potential.acceleration(&x, t);
            for d in 0..3 {
                v[d] += 0.5 * h * a[d];
            }
        }
        // Avoid accumulating round-off in the output times.
        t = t0 + k as f64 * out_dt;
        ts.push(t);
        xs.push(x);
        vs.push(v);
    }
    OrbitTable::new(ts, xs, vs).expect("integrator produces increasing times")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::potentials::{Nfw, PointMass};
    use crate::units::G;
    use std::f64::consts::PI;

    fn energy(p: &dyn Potential, x: &[f64; 3], v: &[f64; 3], t: f64) -> f64 {
        0.5 * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) + p.potential(x, t)
    }

    #[test]
    fn circular_kepler_orbit_closes_after_one_period() {
        let p = PointMass { mass: 1.0.into() };
        let r = 10.0;
        let v_c = (G * 1.0 / r).sqrt();
        let period = 2.0 * PI * r / v_c;
        let orbit = integrate(&p, [r, 0.0, 0.0], [0.0, v_c, 0.0], 0.0, period, period / 2000.0, period / 100.0);
        let (x, v) = orbit.at(period).unwrap();
        assert!(((x[0] - r).powi(2) + x[1].powi(2)).sqrt() < 1e-4 * r, "end position {x:?}");
        let e0 = energy(&p, &[r, 0.0, 0.0], &[0.0, v_c, 0.0], 0.0);
        assert!((energy(&p, &x, &v, period) / e0 - 1.0).abs() < 1e-6);
        assert_eq!(orbit.t.len(), 101);
    }

    #[test]
    fn eccentric_orbit_in_nfw_conserves_energy_and_angular_momentum() {
        let cluster = Nfw::from_virial(1e5, 2000.0, 5.0);
        let x0 = [1000.0, 0.0, 0.0];
        let v0 = [0.0, 300.0, 100.0];
        // Plunging orbit: pericenter about 100 kpc, steps of about 1 Myr.
        let orbit = integrate(&cluster, x0, v0, 0.0, 10.0, 1e-3, 0.05);
        let e0 = energy(&cluster, &x0, &v0, 0.0);
        let l = |x: &[f64; 3], v: &[f64; 3]| [x[1] * v[2] - x[2] * v[1], x[2] * v[0] - x[0] * v[2], x[0] * v[1] - x[1] * v[0]];
        let l0 = l(&x0, &v0);
        let l0n = (l0[0] * l0[0] + l0[1] * l0[1] + l0[2] * l0[2]).sqrt();
        let r_min = orbit.pos.iter().map(|x| (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt()).fold(f64::MAX, f64::min);
        let de_max = (0..orbit.t.len()).map(|i| (energy(&cluster, &orbit.pos[i], &orbit.vel[i], orbit.t[i]) / e0 - 1.0).abs()).fold(0.0, f64::max);
        println!("r_min = {r_min:.1} kpc, max |dE/E| = {de_max:.2e}");
        for i in 0..orbit.t.len() {
            let (x, v) = (orbit.pos[i], orbit.vel[i]);
            assert!((energy(&cluster, &x, &v, orbit.t[i]) / e0 - 1.0).abs() < 1e-4, "energy at row {i}");
            let li = l(&x, &v);
            for d in 0..3 {
                assert!((li[d] - l0[d]).abs() < 1e-9 * l0n, "angular momentum at row {i}");
            }
        }
    }
}
