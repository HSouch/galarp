//! A column-density based shadow model. Fully dynamic and takes into account
//! the disk-wind angle.
//! 

use std::f64::consts::PI;

// use super::utils::{dot, basis, sub, norm};
use super::utils::basis;
use crate::utils::{dot, sub, norm};

use rayon::prelude::*;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DynamicShadow {
    pub bin_width: f64,

    pub refill_length: f64,
    /// Radius where you don't shadow anymore
    pub cutoff_radius: f64,
    pub cutoff_taper: f64,
}


impl Default for DynamicShadow {
    fn default() -> Self {
        Self {
            bin_width: 2.0,
            refill_length: f64::INFINITY,
            cutoff_radius: f64::INFINITY,
            cutoff_taper: 0.0,
        }
    }
}


impl DynamicShadow {
    /// Shadow factoring S of particles at `pos` with cross sections `area`. 
    /// `los` is the unit direction the wind blows along, 
    /// `north` fixes the orientation of the column grid about it (GalaRP uses the disk normal) 
    /// `center` is the origin.
    /// 
    /// Shadowing operates first by binning the particles by x_cam and y_cam, which are
    /// set by the line of sight `los` and rotation `north`. It then computes shadowing on the 
    /// particles in each column by adopting an optical-depth prescription.
    /// 
    pub fn factors(&self, los: [f64; 3], north: [f64; 3], center: [f64; 3], pos: &[[f64; 3]], area: &[f64]) -> Vec<f64> {
        assert_eq!(pos.len(), area.len());
        let (xcam, ycam) = basis(los, north);
        let w = self.bin_width;
        
        // (column x, column y, distance along the wind) per particle.
        let keys: Vec<(i64, i64, f64)> = pos
            .par_iter()
            .map(|x| {
                let d = [x[0] - center[0], x[1] - center[1], x[2] - center[2]];
                ((dot(&d, &xcam) / w).floor() as i64, (dot(&d, &ycam) / w).floor() as i64, dot(&d, &los))
            })
            .collect();
        
        // Ordering of particles within the columns
        let mut order: Vec<usize> = (0..pos.len()).collect();
        order.par_sort_unstable_by(|&a, &b| {
            let (ka, kb) = (&keys[a], &keys[b]);
            (ka.0, ka.1).cmp(&(kb.0, kb.1)).then(ka.2.total_cmp(&kb.2)).then(a.cmp(&b))
        });
        
        let columns: Vec<&[usize]> = 
            order.chunk_by(|&a, &b| keys[a].0 == keys[b].0 && keys[a].1 == keys[b].1).collect();
        
        let weight: Vec<f64> = pos
            .par_iter()
            .map(|x| self.weight(norm(&sub(&x, &center))))
            .collect();
        
        let cell_area = w * w;
        // March through the columns and determine shadowing for those particles.
        // Done in parallel using rayon par_iter()
        let marched: Vec<Vec<(usize, f64)>> = columns
            .par_iter()
            .map(|col| {
                let mut phi = 1.0;
                let mut d_prev: Option<f64> = None;
                col.iter().map(|&i| {
                    let d = keys[i].2;
                    if let Some(dp) = d_prev {
                        phi = 1.0 - (1.0 - phi) * (-(d-dp) / self.refill_length).exp();
                    }
                    d_prev = Some(d);
                    
                    let wt = weight[i];

                    let tau = wt * area[i] / cell_area;

                    // 1 - e^-tau, accurate for small tau.
                    let absorbed = -(-tau).exp_m1();
                    let self_shielding = if tau > 0.0 { absorbed / tau } else { 1.0 };
                    let s = phi * self_shielding;
                    phi *= 1.0 - absorbed;
                    (i, s)
                
                })
                .collect()
            })
            .collect();
        let mut s = vec![1.0; pos.len()];
        for (i, v) in marched.into_iter().flatten() {
            s[i] = v;
        }
        s
    }

    fn weight(&self, r: f64) -> f64 {
        if r <= self.cutoff_radius { 1.0 }
        else if r >= self.cutoff_radius + self.cutoff_taper { 0.0 }
        else { 0.5 * (1.0 + (PI * (r - self.cutoff_radius) / self.cutoff_taper).cos())} 
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn upstream_is_opposite_the_wind_direction() {
        // Wind blowing along -x: the particle at larger x is upstream.
        let shadow = DynamicShadow::default();
        let pos = [[-2.0, 0.2, 0.3], [2.0, 0.2, 0.3]];
        let s = shadow.factors([-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0; 3], &pos, &[0.5; 2]);
        assert!(s[1] > s[0], "{s:?}");
    }

    #[test]
    fn basis_is_orthonormal_even_along_north() {
        for los in [[0.0, 0.0, 1.0], [0.0, 0.0, -1.0], [0.6, 0.0, 0.8], [0.0, 1.0, 0.0]] {
            let (x, y) = basis(los, [0.0, 0.0, 1.0]);
            for (a, b) in [(x, y), (x, los), (y, los)] {
                assert!(dot(&a, &b).abs() < 1e-12);
            }
            assert!((dot(&x, &x) - 1.0).abs() < 1e-12 && (dot(&y, &y) - 1.0).abs() < 1e-12);
        }
    }
}
