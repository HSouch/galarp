//! Handling of external gravitational field acting on the clouds: 
//!     static cluster potential and any number of galaxy potentials
//! 


use crate::orbit::OrbitTable;
use crate::orientation::Orientation;
use crate::potentials::Potential;


pub struct GalaxyField {

    pub potential: Box<dyn Potential>,

    pub orbit: Option<OrbitTable>,

    pub orientation: Orientation,
}


impl GalaxyField {
    pub fn new(potential: impl Potential + 'static, orbit: Option<OrbitTable>) -> Self {
        Self {
            potential: Box::new(potential),
            orbit,
            orientation: Orientation::IDENTITY
        }
    }


    pub fn with_orientation(mut self, orientation: Orientation) -> Self {
        self.orientation = orientation;
        self
    }

    #[inline]
    pub fn acceleration_at_offset(&self, d: &[f64; 3], t: f64) -> [f64; 3] {
        if self.orientation.is_identity() { 
            return self.potential.acceleration(d, t); 
        } else {
            return self.orientation.to_cluster(&self.potential.acceleration(&self.orientation.to_galaxy(d), t));
        }
        
    }

    /// Potential at offset `d` (cluster-frame axes) from the galaxy center.
    #[inline]
    pub fn potential_at_offset(&self, d: &[f64; 3], t: f64) -> f64 {
        self.potential.potential(&self.orientation.to_galaxy(d), t)
    }

    /// Center and velocity in cluster coordinates at time `t`.
    ///
    /// Panics if `t` is outside the orbit table: running past the orbit is a setup error.
    pub fn state(&self, t: f64) -> ([f64; 3], [f64; 3]) {
        match &self.orbit {
            None => ([0.0; 3], [0.0; 3]),
            Some(orbit) => orbit.at(t).unwrap_or_else(|| {
                panic!("time {t} is outside the orbit table [{}, {}]", orbit.t_start(), orbit.t_end())
            }),
        }
    }

}


#[derive(Default)]
pub struct ExternalField {
    pub cluster: Option<Box<dyn Potential>>,

    pub galaxies: Vec<GalaxyField>,
}


impl ExternalField {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_cluster(mut self, p: impl Potential + 'static) -> Self {
        self.cluster = Some(Box::new(p));
        self
    }

    pub fn with_galaxy(mut self, p: impl Potential + 'static, orbit: Option<OrbitTable>) -> Self {
        self.galaxies.push(GalaxyField::new(p, orbit));
        self
    }

    /// Parent galaxy position and velocity in cluster coords at time t
    pub fn galaxy_state(&self, t: f64) -> ([f64; 3], [f64; 3]) {
        self.galaxies.first().map_or(([0.0; 3], [0.0; 3]), |g| g.state(t))
    }

    /// Centers of all galaxies at time `t`, in cluster coordinates.
    pub fn galaxy_centers(&self, t: f64) -> Vec<[f64; 3]> {
        self.galaxies.iter().map(|g| g.state(t).0).collect()
    }

    /// Acceleration at cluster-frame position `xc` and time `t`, given the galaxy centers at
    /// that time (from `galaxy_centers`, computed once per step rather than per cloud).
    #[inline]
    pub fn acceleration_with_centers(&self, xc: &[f64; 3], t: f64, centers: &[[f64; 3]]) -> [f64; 3] {
        let mut a = [0.0; 3];
        if let Some(p) = &self.cluster {
            a = p.acceleration(xc, t);
        }
        for (g, c) in self.galaxies.iter().zip(centers) {
            let ag = g.acceleration_at_offset(&[xc[0] - c[0], xc[1] - c[1], xc[2] - c[2]], t);
            for k in 0..3 {
                a[k] += ag[k];
            }
        }
        a
    }

    pub fn acceleration(&self, x: &[f64; 3], t: f64) -> [f64; 3] {
        self.acceleration_with_centers(x, t, &self.galaxy_centers(t))
    }

}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::potentials::{Hernquist, Nfw, PointMass};


    #[test]
    fn several_galaxies_add_and_the_first_is_the_parent() {
        let line = |x0: f64, v: f64| {
            OrbitTable::new(vec![0.0, 1.0], vec![[x0, 0.0, 0.0], [x0 + v, 0.0, 0.0]], vec![[v, 0.0, 0.0], [v, 0.0, 0.0]]).unwrap()
        };

        let a = PointMass {mass: 1.0.into() };
        let b = PointMass { mass: 2.0.into() };
        let field = ExternalField::new()
            .with_galaxy(a.clone(), Some(line(0.0, 10.0)))
            .with_galaxy(b.clone(), Some(line(100.0, -10.0)));
        assert_eq!(field.galaxy_state(0.5).0, [5.0, 0.0, 0.0]);

        let x = [50.0, 3.0, 0.0];
        let sum = {
            let aa = a.acceleration(&[45.0, 3.0, 0.0], 0.5);
            let bb = b.acceleration(&[-45.0, 3.0, 0.0], 0.5);
            [aa[0] + bb[0], aa[1] + bb[1], aa[2] + bb[2]]
        };
        assert_eq!(field.acceleration(&x, 0.5), sum);
    }

    #[test]
    fn galaxy_potential_follows_the_orbit_and_cluster_sits_at_the_origin() {
        let orbit = OrbitTable::new(
            vec![0.0, 1.0],                                 // Times
            vec![[100.0, 0.0, 0.0], [200.0, 0.0, 0.0]],     // Positions
            vec![[100.0, 0.0, 0.0], [100.0, 0.0, 0.0]],     // Velocities
        )
        .unwrap();
        
        let galaxy = Hernquist { mass: 10.0.into(), a: 2.0.into() };
        let field = ExternalField::new().with_galaxy(galaxy.clone(), Some(orbit));

        // At t = 0.5 the galaxy is at x = 150.
        let at_center = field.acceleration(&[150.0, 0.0, 0.0], 0.5);
        assert!(at_center.iter().all(|v| v.abs() < 1e-12));
        let a = field.acceleration(&[153.0, 0.0, 4.0], 0.5);
        let expected = galaxy.acceleration(&[3.0, 0.0, 4.0], 0.5);
        assert_eq!(a, expected);

        let cluster = Nfw { m: 1e4.into(), r_s: 300.0.into() };
        let field = ExternalField::new().with_cluster(cluster.clone());
        assert_eq!(field.acceleration(&[510.0, 20.0, 30.0], 0.0), cluster.acceleration(&[510.0, 20.0, 30.0], 0.0));
    }
}