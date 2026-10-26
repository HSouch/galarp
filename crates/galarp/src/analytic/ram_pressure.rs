

use rayon::prelude::*;

use super::utils::{dot, norm};
use super::shadow::DynamicShadow;

use crate::cluster::ClusterGas;
use crate::external::ExternalField;
use crate::potentials::Param;


/// ICM density calculated at position
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum DensityAt {
    #[default]
    Galaxy,
    Particle,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Wind {
    pub center: [f64; 3],
    pub velocity: [f64; 3],

    pub normal: [f64; 3],
    pub density: f64,
}


impl Wind {
    pub fn speed(&self) -> f64 {
        norm(&self.velocity)
    }

    pub fn direction(&self) -> Option<[f64; 3]> {
        let s = self.speed();
        (s > 0.0).then(|| [-self.velocity[0] / s, -self.velocity[1] / s, -self.velocity[2] / s])
    }

    pub fn ram_pressure(&self) -> f64 {
        self.density * self.speed().powi(2)
    }

    pub fn angle_deg(&self) -> f64 {
        match self.direction() {
            Some(w) => dot(&w, &self.normal).abs().min(1.0).acos().to_degrees(),
            None => f64::NAN,
        }
    }
}


pub struct RamPressure {
    pub icm: ClusterGas,
    pub surface_density: Param,
    pub density_at: DensityAt,
    pub shadow: Option<DynamicShadow>,
}


impl RamPressure {
    pub fn winds(&self, field: &ExternalField, t: f64) -> Vec<Wind> {
        field
            .galaxies
            .iter()
            .map(|g| {
                let (center, velocity) = g.state(t);
                Wind { 
                        center, 
                        velocity, 
                        normal: g.orientation.normal(), 
                        density: self.icm.density_at(&center)
                    }
            })
            .collect()
    }

    pub fn shadowing(&self, winds: &[Wind],  t: f64, pos: &[[f64; 3]], mass: &[f64], 
        galaxy: &[u16]) -> Vec<f64> {
            let mut s = vec![1.0; pos.len()];
            
            // If no shadow implemented, just return 1.0 for all clouds
            let Some(shadow) = &self.shadow else { return s };

            let sigma = self.surface_density.at(t);

            for (k, wind) in winds.iter().enumerate() {
                let Some(los) = wind.direction() else { continue };
                let members: Vec<usize> = (0..pos.len())
                    .filter(|&i| galaxy[i] as usize == k).collect();
                    if members.is_empty() { continue; }
                let x: Vec<[f64; 3]> = members.iter()
                    .map(|&i| pos[i]).collect();
                let area: Vec<f64> = members.iter()
                    .map(|&i| mass[i] / sigma).collect();
                
                for (&i, f) in members.iter()
                    .zip(shadow.factors(los, wind.normal, wind.center, &x, &area)) {
                        s[i] = f;
                    }
                }
            s
        }
    
    /// Gets acceleration for a given cloud
    pub fn acceleration(&self, wind: &Wind, sigma: f64, x: &[f64; 3], v: &[f64; 3], s: f64) 
        -> [f64; 3] {
            let rho = match self.density_at {
                DensityAt::Galaxy => wind.density,
                DensityAt::Particle => self.icm.density_at(x),
            };

            let v_rel = [-v[0], -v[1], -v[2]];
            let f = s * rho * norm(&v_rel) / sigma;

            [f * v_rel[0], f * v_rel[1], f * v_rel[2]]
        }


    /// Computes acceleration and shadowing for every cloud particle
    pub fn evaluate(&self, field: &ExternalField, t: f64, pos: &[[f64; 3]], vel: &[[f64; 3]],
        mass: &[f64], galaxy: &[u16]) -> (Vec<[f64; 3]>, Vec<f64>) {
            let winds = self.winds(field, t);

            let s = self.shadowing(&winds, t, pos, mass, galaxy);

            let sigma = self.surface_density.at(t);
            let acc = (0..pos.len())
                .into_par_iter()
                .map(|i| self.acceleration(
                    &winds[galaxy[i] as usize],
                    sigma, 
                    &pos[i],
                    &vel[i],
                    s[i]))
                .collect();
           
            (acc, s)
            
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::orbit::OrbitTable;
    use crate::potentials::PointMass;

    fn moving_galaxy(v: [f64; 3]) -> ExternalField {
        let orbit = OrbitTable::new(vec![0.0, 1.0], vec![[0.0; 3], v], vec![v, v]).unwrap();
        ExternalField::new().with_galaxy(PointMass { mass: 0.0.into() }, Some(orbit))
    }

    #[test]
    fn acceleration_is_rho_v_squared_over_sigma_against_the_motion() {
        let rp = RamPressure { icm: ClusterGas::Uniform { rho: 2.0, u: 1.0 }, surface_density: 4.0.into(), density_at: DensityAt::Galaxy, shadow: None };
        let field = moving_galaxy([0.0, 0.0, 3.0]);
        let wind = rp.winds(&field, 0.5)[0];
        assert_eq!(wind.direction(), Some([0.0, 0.0, -1.0]));
        assert_eq!(wind.ram_pressure(), 18.0);
        // A cloud moving with the galaxy plus 4 sideways: |v| = 5.
        let a = rp.acceleration(&wind, 4.0, &[0.0; 3], &[4.0, 0.0, 3.0], 1.0);
        let expected = [-2.0 * 5.0 * 4.0 / 4.0, 0.0, -2.0 * 5.0 * 3.0 / 4.0];
        assert_eq!(a, expected);
        assert_eq!(wind.angle_deg(), 0.0);
    }

    #[test]
    fn each_galaxy_shadows_only_its_own_clouds() {
        let rp = RamPressure {
            icm: ClusterGas::Uniform { rho: 1.0, u: 1.0 },
            surface_density: 1.0.into(),
            density_at: DensityAt::Galaxy,
            shadow: Some(DynamicShadow { drag_coeff: 1.0, bin_width: 1.0 }),
        };
        // Two galaxies, both at z = 0.5 moving along +z: the wind blows along -z, so the
        // cloud above is upstream. The third cloud is in the same column but belongs to the
        // second galaxy, so it is alone in its own galaxy's columns and unshadowed.
        let orbit = OrbitTable::new(vec![0.0, 1.0], vec![[0.0; 3], [0.0, 0.0, 1.0]], vec![[0.0, 0.0, 1.0]; 2]).unwrap();
        let field = moving_galaxy([0.0, 0.0, 1.0]).with_galaxy(PointMass { mass: 0.0.into() }, Some(orbit));
        let pos = [[0.2, 0.2, 2.0], [0.2, 0.2, 1.0], [0.2, 0.2, 0.0]];
        let (acc, s) = rp.evaluate(&field, 0.5, &pos, &[[0.0, 0.0, 1.0]; 3], &[0.5; 3], &[0, 0, 1]);
        assert!(s[0] < 1.0 && s[1] < s[0] && s[2] == 1.0, "{s:?}");
        // Each cloud moves with its galaxy, so it meets a unit wind along -z scaled by S.
        assert!(acc.iter().zip(&s).all(|(a, &si)| (a[2] + si).abs() < 1e-12), "{acc:?}");
    }
}
