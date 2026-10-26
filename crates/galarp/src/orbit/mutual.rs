//! Integrates several galaxies together, each pulled by the cluster potential and by every
//! other galaxy's potential.
//!


use super::table::OrbitTable;
use crate::orientation::Orientation;
use crate::potentials::Potential;

/// How one galaxy moves.
pub enum Motion {
    /// Integrate from position and velocity at t0 (cluster frame).
    Integrate { x0: [f64; 3], v0: [f64; 3] },
    /// Follow this table (it must cover [t0, t_end]).
    Prescribed(OrbitTable),
}

/// One galaxy: its potential, its orientation (galaxy frame to cluster frame) and motion.
pub struct Body<'a> {
    pub potential: &'a dyn Potential,
    pub orientation: Orientation,
    pub motion: Motion,
}

impl Body<'_> {
    /// Acceleration this galaxy's potential gives at offset `d` (cluster axes) from its center.
    fn pull(&self, d: &[f64; 3], t: f64) -> [f64; 3] {
        if self.orientation.is_identity() {
            return self.potential.acceleration(d, t);
        }
        self.orientation.to_cluster(&self.potential.acceleration(&self.orientation.to_galaxy(d), t))
    }
}

/// Acceleration of every integrated body at time `t`, with the bodies at `x` (prescribed
/// ones already placed at their table positions).
fn accelerations(cluster: Option<&dyn Potential>, bodies: &[Body], x: &[[f64; 3]], mutual: bool, t: f64) -> Vec<[f64; 3]> {
    (0..bodies.len())
        .map(|j| {
            if matches!(bodies[j].motion, Motion::Prescribed(_)) {
                return [0.0; 3];
            }
            let mut a = cluster.map_or([0.0; 3], |c| c.acceleration(&x[j], t));
            if mutual {
                for (i, other) in bodies.iter().enumerate() {
                    if i == j {
                        continue;
                    }
                    let d = [x[j][0] - x[i][0], x[j][1] - x[i][1], x[j][2] - x[i][2]];
                    let ai = other.pull(&d, t);
                    for k in 0..3 {
                        a[k] += ai[k];
                    }
                }
            }
            a
        })
        .collect()
}

/// Integrates the bodies from `t0` to `t_end`, returning one orbit table per body (a
/// prescribed body's own table, unchanged). Rows every `output_dt` (rounded to a whole
/// number of rows), each split into equal steps no longer than `dt`. With `mutual` false,
/// each integrated body feels only the cluster, as with `integrate`.
pub fn integrate_bodies(
    cluster: Option<&dyn Potential>,
    bodies: &[Body],
    mutual: bool,
    t0: f64,
    t_end: f64,
    dt: f64,
    output_dt: f64,
) -> Vec<OrbitTable> {
    assert!(t_end > t0 && dt > 0.0 && output_dt > 0.0, "need t_end > t0 and positive steps");
    let n_out = ((t_end - t0) / output_dt).round().max(1.0) as usize;
    let out_dt = (t_end - t0) / n_out as f64;
    let n_sub = (out_dt / dt).ceil().max(1.0) as usize;
    let h = out_dt / n_sub as f64;

    let place = |b: &Body, t: f64| match &b.motion {
        Motion::Integrate { .. } => None,
        Motion::Prescribed(table) => Some(table.at(t).unwrap_or_else(|| panic!("prescribed orbit does not cover t = {t}"))),
    };
    let mut x: Vec<[f64; 3]> = bodies.iter().map(|b| match &b.motion {
        Motion::Integrate { x0, .. } => *x0,
        Motion::Prescribed(table) => table.at(t0).expect("prescribed orbit covers t0").0,
    }).collect();
    let mut v: Vec<[f64; 3]> = bodies.iter().map(|b| match &b.motion {
        Motion::Integrate { v0, .. } => *v0,
        Motion::Prescribed(table) => table.at(t0).expect("prescribed orbit covers t0").1,
    }).collect();

    let mut t = t0;
    let mut a = accelerations(cluster, bodies, &x, mutual, t);
    let mut ts = vec![t];
    let mut xs: Vec<Vec<[f64; 3]>> = x.iter().map(|&p| vec![p]).collect();
    let mut vs: Vec<Vec<[f64; 3]>> = v.iter().map(|&p| vec![p]).collect();
    for k in 1..=n_out {
        for _ in 0..n_sub {
            for j in 0..bodies.len() {
                if matches!(bodies[j].motion, Motion::Integrate { .. }) {
                    for d in 0..3 {
                        v[j][d] += 0.5 * h * a[j][d];
                        x[j][d] += h * v[j][d];
                    }
                }
            }
            t += h;
            for (j, b) in bodies.iter().enumerate() {
                if let Some((xp, vp)) = place(b, t) {
                    (x[j], v[j]) = (xp, vp);
                }
            }
            a = accelerations(cluster, bodies, &x, mutual, t);
            for j in 0..bodies.len() {
                if matches!(bodies[j].motion, Motion::Integrate { .. }) {
                    for d in 0..3 {
                        v[j][d] += 0.5 * h * a[j][d];
                    }
                }
            }
        }
        // Avoid accumulating round-off in the output times.
        t = t0 + k as f64 * out_dt;
        ts.push(t);
        for j in 0..bodies.len() {
            xs[j].push(x[j]);
            vs[j].push(v[j]);
        }
    }

    bodies
        .iter()
        .zip(xs.into_iter().zip(vs))
        .map(|(b, (xj, vj))| match &b.motion {
            Motion::Prescribed(table) => table.clone(),
            Motion::Integrate { .. } => OrbitTable::new(ts.clone(), xj, vj).expect("integrator produces increasing times"),
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::orbit::integrate;
    use crate::potentials::{Hernquist, Nfw, PointMass};
    use crate::units::G;
    use std::f64::consts::PI;

    fn body(p: &dyn Potential, x0: [f64; 3], v0: [f64; 3]) -> Body<'_> {
        Body { potential: p, orientation: Orientation::IDENTITY, motion: Motion::Integrate { x0, v0 } }
    }

    #[test]
    fn one_body_matches_the_single_orbit_integrator_exactly() {
        let cluster = Nfw::from_virial(1e5, 2000.0, 5.0);
        let g = Hernquist { mass: 1.0.into(), a: 1.0.into() };
        let (x0, v0) = ([1000.0, 0.0, 0.0], [0.0, 300.0, 100.0]);
        let single = integrate(&cluster, x0, v0, 0.0, 5.0, 1e-3, 0.05);
        let both = integrate_bodies(Some(&cluster), &[body(&g, x0, v0)], true, 0.0, 5.0, 1e-3, 0.05);
        assert_eq!(single.pos, both[0].pos);
        assert_eq!(single.vel, both[0].vel);
    }

    #[test]
    fn two_point_masses_orbit_their_center_of_mass() {
        // Equal point masses on a circular binary: separation 2r, each moving at
        // v = sqrt(G m / (4 r)) about the fixed center of mass.
        let m = 1.0;
        let r = 10.0;
        let v = (G * m / (4.0 * r)).sqrt();
        let p = PointMass { mass: m.into() };
        let period = 2.0 * PI * r / v;
        let bodies = [body(&p, [r, 0.0, 0.0], [0.0, v, 0.0]), body(&p, [-r, 0.0, 0.0], [0.0, -v, 0.0])];
        let orbits = integrate_bodies(None, &bodies, true, 0.0, period, period / 4000.0, period / 50.0);
        for i in 0..orbits[0].t.len() {
            let (a, b) = (orbits[0].pos[i], orbits[1].pos[i]);
            // The center of mass stays put and the separation stays 2r.
            assert!((0..3).all(|k| (a[k] + b[k]).abs() < 1e-9 * r), "row {i}");
            let sep = ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt();
            assert!((sep / (2.0 * r) - 1.0).abs() < 1e-5, "row {i}: separation {sep}");
        }
        // Back where it started after one period.
        let end = *orbits[0].pos.last().unwrap();
        assert!((end[0] - r).abs() < 1e-3 * r && end[1].abs() < 1e-3 * r, "{end:?}");
        // Without mutual gravity (and no cluster) they fly apart in straight lines.
        let free = integrate_bodies(None, &bodies, false, 0.0, period, period / 4000.0, period / 50.0);
        let end = free[0].pos.last().unwrap();
        assert!((end[1] - v * period).abs() < 1e-9 * v * period);
    }

    #[test]
    fn a_prescribed_body_pulls_but_is_not_pulled() {
        // A satellite on a circular orbit about a host held fixed at the origin by its table.
        let host_p = PointMass { mass: 1.0.into() };
        let sat_p = PointMass { mass: 0.5.into() };
        let r = 20.0;
        let v = (G * 1.0 / r).sqrt();
        let period = 2.0 * PI * r / v;
        let fixed = OrbitTable::new(vec![0.0, 2.0 * period], vec![[0.0; 3]; 2], vec![[0.0; 3]; 2]).unwrap();
        let bodies = [
            Body { potential: &host_p, orientation: Orientation::IDENTITY, motion: Motion::Prescribed(fixed.clone()) },
            body(&sat_p, [r, 0.0, 0.0], [0.0, v, 0.0]),
        ];
        let orbits = integrate_bodies(None, &bodies, true, 0.0, period, period / 4000.0, period / 40.0);
        assert_eq!(orbits[0].pos, fixed.pos, "prescribed table returned unchanged");
        for x in &orbits[1].pos {
            let rr = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
            assert!((rr / r - 1.0).abs() < 1e-5, "radius {rr}");
        }
    }
}
