

use std::f64::consts::PI;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum GridLayout { Square, Rings }

#[derive(Clone, Copy, Debug)]
pub struct TestGrid {
    pub r_max: f64,
    pub spacing: f64,
    /// Gas surface density, in code mass per kpc^2.
    pub sigma_gas: f64,
    pub layout: GridLayout,
}


impl TestGrid {
    /// Positions (galaxy frame, z = 0) and the common cloud mass.
    pub fn points(&self) -> (Vec<[f64; 3]>, f64) {
        let d = self.spacing;
        let mut pos = Vec::new();
        match self.layout {
            GridLayout::Square => {
                let n = (self.r_max / d).ceil() as i64;
                for i in -n..n {
                    for j in -n..n {
                        let (x, y) = ((i as f64 + 0.5) * d, (j as f64 + 0.5) * d);
                        
                        // Only add the position if it's within r_max
                        if x.hypot(y) < self.r_max {
                            pos.push([x, y, 0.0]);
                        }
                    }
                }
            }
            GridLayout::Rings => {
                let n_rings = (self.r_max / d).floor() as usize;
                for k in 0..n_rings {
                    let r = (k as f64 + 0.5) * d;
                    let n_k = ((2.0 * PI * r / d).round() as usize).max(1);
                    let offset = if k % 2 == 1 { 0.5 } else { 0.0 };
                    for j in 0..n_k {
                        let phi = 2.0 * PI * (j as f64 + offset) / n_k as f64;
                        pos.push([r * phi.cos(), r * phi.sin(), 0.0]);
                    }
                }
            }
        }
        (pos, self.sigma_gas * d * d)
    }
}
