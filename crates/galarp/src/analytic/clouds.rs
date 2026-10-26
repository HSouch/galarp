//! Cloud test particles and ICs

use crate::ic::GasDisk;

use crate::orientation::Orientation;
use crate::potentials::Potential;
use crate::rng::Rng;


/// Cloud particles
/// Positions and velocities are in cluster coordinates.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Clouds {
    pub id: Vec<u64>,
    pub galaxy: Vec<u16>,
    pub pos: Vec<[f64; 3]>,
    pub vel: Vec<[f64; 3]>,
    pub mass: Vec<f64>,

    pub shadow: Vec<f64>,
    pub a_ram: Vec<[f64; 3]>,
}


impl Clouds {
    pub fn len(&self) -> usize {
        self.id.len()
    }

    pub fn is_empty(&self) -> bool {
        self.id.is_empty()
    }

    pub fn push(&mut self, id: u64, galaxy: u16, pos: [f64; 3], vel: [f64; 3], mass: f64) {
        self.id.push(id);
        self.galaxy.push(galaxy);
        self.pos.push(pos);
        self.vel.push(vel);
        self.mass.push(mass);
        self.shadow.push(1.0);
        self.a_ram.push([0.0; 3])
    }

    pub fn add_disk(
        &mut self,
        disk: &GasDisk,
        potential: &dyn Potential,
        orientation: Orientation,
        galaxy: u16,
        n: usize,
        center: [f64; 3],
        bulk_velocity: [f64; 3],
        sigma_z: f64,
        seed: u64,
        next_id: &mut u64,
    ) {
        
        // Individual cloud masses defined by disk total mass and number of clouds.
        let m = disk.params.mass / n as f64;
        let mut rng = Rng::new(seed);

        // Sample the disk n times.
        for x in disk.sample(n) {
            let r = x[0].hypot(x[1]);

            let v_phi = 
                if r > 0.0 {
                    let a = potential.acceleration(&x, 0.0);
                    (-(a[0] * x[0] + a[1] * x[1])).max(0.0).sqrt()
                } else {
                    0.0
                };
            
            let (s, c) = 
                if r > 0.0 { (x[1] / r, x[0] / r) } 
                else { (0.0, 1.0) };
        
            
            let vz: f64 = sigma_z * rng.sample::<f64, _>(&rand_distr::StandardNormal);

            let (xc, vc) = (
                orientation.to_cluster(&x), 
                orientation.to_cluster(&[-v_phi * s,  v_phi * c, vz])
            );

            self.push(
                *next_id,
                galaxy,
                [xc[0] + center[0], xc[1] + center[1], xc[2] + center[2]],
                [vc[0] + bulk_velocity[0], vc[1] + bulk_velocity[1], vc[2] + bulk_velocity[2]],
                m,
            );
            *next_id += 1;
        }

    }
}


#[cfg(test)]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::ic::GasDiskParams;
    use crate::potentials::{Composite, Hernquist, MiyamotoNagai};
    use crate::units::G;

    #[test]
    fn disk_clouds_rotate_at_the_circular_velocity_and_carry_the_disk_mass() {
        let pot = Composite::new()
            .with(MiyamotoNagai { mass: 5.0.into(), a: 3.0.into(), b: 0.3.into() })
            .with(Hernquist { mass: 1.0.into(), a: 0.5.into() });
        
        let disk = GasDisk::new(GasDiskParams { mass: 0.5, scale_length: 4.0, r_max: 20.0, sound_speed: 20.0 }, &pot);
        let mut clouds = Clouds::default();
        let mut id = 0;
        let o = Orientation::from_angles(30.0, 0.0);
        
        clouds.add_disk(&disk, &pot, o, 0, 500, [100.0, 0.0, 0.0], [0.0, 0.0, 50.0], 0.0, 1, &mut id);
        assert_eq!(clouds.len(), 500);
        assert!((clouds.mass.iter().sum::<f64>() - 0.5).abs() < 1e-12);
        
        for i in 0..clouds.len() {
            let d = [clouds.pos[i][0] - 100.0, clouds.pos[i][1], clouds.pos[i][2]];
            let v = [clouds.vel[i][0], clouds.vel[i][1], clouds.vel[i][2] - 50.0];
            let (xg, vg) = (o.to_galaxy(&d), o.to_galaxy(&v));
            let r = xg[0].hypot(xg[1]);
            // In the disk plane, perpendicular to the radius, and circular at large radius.
            assert!(vg[2].abs() < 1e-9 && (xg[0] * vg[0] + xg[1] * vg[1]).abs() < 1e-9 * r * 300.0);
            if r > 10.0 && xg[2].abs() < 0.05 {
                let v_kepler = (G * 6.0 / r).sqrt();
                assert!((vg[0].hypot(vg[1]) / v_kepler - 1.0).abs() < 0.15);
            }
        }
    }
}
