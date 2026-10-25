


//! Isotherman gas disk in HSE in a static potential
//! 
//! Optically can include the gas disk's own gravity.
//! 
//! Surface density: Sigma(R) = Sigma_0 exp(-R / R_d) for R < R_max
//! 
//! Vertical Structure: At each radius, HSE Gives:
//! 
//! ```text
//! rho(R, z) = rho_0(R) exp(-[dPhi_ext(R, z) + dPhi_gas(R, z)] / s_s^2)
//! ```
//!
//! To avoid clumping, we use Halton sampling.
//! 


use std::f64::consts::PI;

use rayon::prelude::*;

use super::halton;
use crate::potentials::Potential;
use crate::rng::Rng;
use crate::units::G;

#[derive(Clone, Copy, Debug)]
pub struct GasDiskParams {
    /// Total gas mass
    pub mass: f64,
    /// Exponential scale length 
    pub scale_length: f64,
    /// Outer truncation radius.
    pub r_max: f64,
    /// Isothermal sound speed 
    pub sound_speed: f64,
}

/// Optional refinements of disk initial conditions.
#[derive(Clone, Copy, Debug, Default)]
pub struct DiskOptions {
    /// Include the gas disk's own gravity in vertical structure
    pub self_gravity: bool,
    /// Gaussian position noise
    pub noise_fraction: f64,
    /// Seed for the noise.
    pub noise_seed: u64,
}

const N_R: usize = 512;
const N_Z: usize = 400;
/// Vertical extent
const Z_DEPTH: f64 = 30.0;

/// Vertical profile at one radius on a uniform grid from 0 to `top`.
struct Row {
    top: f64,
    /// Density shape, normalized so that its integral over all z is 1.
    shape: Vec<f64>,
    /// Cumulative integral of `shape` from 0 (reaches 0.5 at `top`).
    cum: Vec<f64>,
}

impl Row {
    fn dz(&self) -> f64 {
        self.top / (N_Z - 1) as f64
    }

    fn shape_at(&self, z: f64) -> f64 {
        let z = z.abs();
        if z >= self.top {
            return 0.0;
        }
        let f = z / self.dz();
        let k = (f as usize).min(N_Z - 2);
        let w = f - k as f64;
        self.shape[k] * (1.0 - w) + self.shape[k + 1] * w
    }

    /// Height below which a fraction `q` of the upper half's mass lies.
    fn quantile(&self, q: f64) -> f64 {
        let target = q * self.cum[N_Z - 1];
        let j = self.cum
            .partition_point(|&c| c < target)
            .clamp(1, N_Z - 1);
        let w = (target - self.cum[j - 1]) / (self.cum[j] - self.cum[j - 1]).max(f64::MIN_POSITIVE);
        (j as f64 - 1.0 + w) * self.dz()
    }
}

pub struct GasDisk<'a> {
    pub params: GasDiskParams,
    pub options: DiskOptions,
    potential: &'a dyn Potential,
    sigma_0: f64,
    /// Vertical profiles at radii r_k = r_max k / (N_R - 1).
    rows: Vec<Row>,
}

impl<'a> GasDisk<'a> {
    /// A disk in the external potential only, without noise.
    pub fn new(params: GasDiskParams, potential: &'a dyn Potential) -> Self {
        Self::with_options(params, potential, DiskOptions::default())
    }

    pub fn with_options(params: GasDiskParams, potential: &'a dyn Potential, options: DiskOptions) -> Self {
        let x_max = params.r_max / params.scale_length;
        let sigma_0 = params.mass / (2.0 * PI * params.scale_length.powi(2) * (1.0 - (1.0 + x_max) * (-x_max).exp()));
        let mut disk = GasDisk { params, options, potential, sigma_0, rows: Vec::new() };
        let radii: Vec<f64> = (0..N_R).map(|i| params.r_max * i as f64 / (N_R - 1) as f64).collect();
        disk.rows = radii
            .par_iter()
            .map(|&r| disk.vertical_row(r))
            .collect();
        disk
    }

    pub fn sigma(&self, r: f64) -> f64 {
        if r > self.params.r_max {
            0.0
        } else {
            self.sigma_0 * (-r / self.params.scale_length).exp()
        }
    }

    /// Phi_ext(R, z) - Phi_ext(R, 0), over c^2.
    fn ext_exponent(&self, r: f64, z: f64) -> f64 {
        let c2 = self.params.sound_speed * self.params.sound_speed;
        (self.potential.potential(&[r, 0.0, z], 0.0) - self.potential.potential(&[r, 0.0, 0.0], 0.0)) / c2
    }

    /// Height where the external Boltzmann factor has dropped to exp(-Z_DEPTH). The gas's
    /// own gravity only makes the disk thinner, so this bounds both cases.
    fn z_top(&self, r: f64) -> f64 {
        let mut z = 1e-4 * self.params.scale_length;
        while self
        .ext_exponent(r, z) < Z_DEPTH && z < 100.0 * self.params.r_max {
            z *= 1.5;
        }
        z
    }

    fn vertical_row(&self, r: f64) -> Row {
        let top = self.z_top(r);
        let dz = top / (N_Z - 1) as f64;
        let e_ext: Vec<f64> = (0..N_Z).map(|k| self.ext_exponent(r, k as f64 * dz)).collect();
        let integrate = |f: &[f64]| {
            let mut cum = vec![0.0; N_Z];
            for k in 1..N_Z {
                cum[k] = cum[k - 1] + 0.5 * (f[k] + f[k - 1]) * dz;
            }
            cum
        };
        let mut f: Vec<f64> = e_ext.iter().map(|e| (-e).exp()).collect();
        if self.options.self_gravity && r < self.params.r_max {
            let sigma = self.sigma(r);
            let c2 = self.params.sound_speed * self.params.sound_speed;
            for _ in 0..200 {
                let cum = integrate(&f);
                let norm = 2.0 * cum[N_Z - 1];

                // Column within +-z, then the slab potential difference dPhi = int 2 pi G Sigma(<z) dz.
                let sigma_in: Vec<f64> = cum.iter().map(|c| sigma * 2.0 * c / norm).collect();
                let dphi = integrate(&sigma_in.iter().map(|s| 2.0 * PI * G * s).collect::<Vec<_>>());
                let new: Vec<f64> = (0..N_Z).map(|k| (-(e_ext[k] + dphi[k] / c2)).exp()).collect();
                let change = (0..N_Z).map(|k| (new[k] - f[k]).abs()).fold(0.0, f64::max);
                
                for k in 0..N_Z {
                    f[k] = 0.5 * (f[k] + new[k]);
                }
                if change < 1e-8 {
                    break;
                }
            }
        }
        let cum = integrate(&f);
        let norm = 2.0 * cum[N_Z - 1];
        Row { top, shape: f.iter().map(|v| v / norm).collect(), cum: cum.iter().map(|c| c / norm).collect() }
    }

    /// Gas density at a position in the galaxy frame.
    pub fn density(&self, x: &[f64; 3]) -> f64 {
        let r = x[0].hypot(x[1]);
        if r >= self.params.r_max {
            return 0.0;
        }
        let f = r / self.params.r_max * (N_R - 1) as f64;
        let k = (f as usize).min(N_R - 2);
        let w = f - k as f64;
        self.sigma(r) * (self.rows[k].shape_at(x[2]) * (1.0 - w) + self.rows[k + 1].shape_at(x[2]) * w)
    }

    /// Mean |z| of the gas at radius `r`.
    pub fn mean_abs_z(&self, r: f64) -> f64 {
        let f = (r / self.params.r_max * (N_R - 1) as f64).min((N_R - 1) as f64);
        let row = &self.rows[f.round() as usize];
        let dz = row.dz();
        let num: f64 = (1..N_Z).map(|k| 0.5 * ((k as f64) * row.shape[k] + (k as f64 - 1.0) * row.shape[k - 1]) * dz * dz).sum();
        num / row.cum[N_Z - 1]
    }

    /// Radius enclosing a fraction `q` of the disk mass.
    fn radius_at_mass_fraction(&self, q: f64) -> f64 {
        let rd = self.params.scale_length;
        let x_max = self.params.r_max / rd;
        let target = q * (1.0 - (1.0 + x_max) * (-x_max).exp());
        let (mut lo, mut hi) = (0.0, x_max);
        for _ in 0..80 {
            let mid = 0.5 * (lo + hi);
            if 1.0 - (1.0 + mid) * (-mid).exp() < target {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi) * rd
    }

    /// Positions (galaxy frame) of `n` equal-mass points, drawn with a Halton sequence
    /// (base 2: enclosed-mass fraction, 3: azimuth, 5: height), then displaced by the
    /// optional noise.
    pub fn sample(&self, n: usize) -> Vec<[f64; 3]> {
        let mut pos: Vec<[f64; 3]> = (0..n)
            .into_par_iter()
            .map(|i| {
                let k = i as u64 + 1;
                let r = self.radius_at_mass_fraction(halton(k, 2));
                let phi = 2.0 * PI * halton(k, 3);
                let uz = halton(k, 5);
                // Fold the vertical sample: lower half of uz below the midplane.
                let (sign, q) = if uz < 0.5 { (-1.0, 1.0 - 2.0 * uz) } else { (1.0, 2.0 * uz - 1.0) };
                let f = (r / self.params.r_max * (N_R - 1) as f64).min((N_R - 1) as f64 - 1e-9);
                let j = (f as usize).min(N_R - 2);
                let w = f - j as f64;
                let z = sign * (self.rows[j].quantile(q) * (1.0 - w) + self.rows[j + 1].quantile(q) * w);
                let (s, c) = phi.sin_cos();
                [r * c, r * s, z]
            })
            .collect();
        if self.options.noise_fraction > 0.0 {
            let m = self.params.mass / n as f64;
            let mut rng = Rng::new(self.options.noise_seed);
            let mut gauss = || -> f64 { rng.sample(&rand_distr::StandardNormal) };
            for x in pos.iter_mut() {
                let spacing = (m / self.density(x).max(f64::MIN_POSITIVE)).cbrt();
                let sigma = self.options.noise_fraction * spacing;
                for c in x.iter_mut() {
                    *c += sigma * gauss();
                }
            }
        }
        pos
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::potentials::{Composite, Hernquist, MiyamotoNagai, Nfw};

    /// Milky Way-like potential: 5e10 Msun disk, 1e10 bulge, NFW halo.
    pub fn milky_way() -> Composite {
        Composite::new()
            .with(MiyamotoNagai { mass: 5.0.into(), a: 3.0.into(), b: 0.3.into() })
            .with(Hernquist { mass: 1.0.into(), a: 0.5.into() })
            .with(Nfw { m: 60.0.into(), r_s: 15.0.into() })
    }

    fn params() -> GasDiskParams {
        GasDiskParams { mass: 0.5, scale_length: 4.0, r_max: 20.0, sound_speed: 20.0 }
    }

    #[test]
    fn density_integrates_to_the_surface_density_and_mass() {
        let mw = milky_way();
        let disk = GasDisk::new(params(), &mw);
        for r in [1.0, 5.0, 12.0] {
            let top = disk.z_top(r);
            let nz = 4000;
            let dz = 2.0 * top / nz as f64;
            let sigma: f64 = (0..nz).map(|k| disk.density(&[r, 0.0, -top + (k as f64 + 0.5) * dz]) * dz).sum();
            assert!((sigma / disk.sigma(r) - 1.0).abs() < 2e-3, "Sigma at R = {r}");
        }
    }

    #[test]
    fn samples_follow_the_surface_density_and_vertical_profile() {
        let mw = milky_way();
        let disk = GasDisk::new(params(), &mw);
        let n = 50_000;
        let parts = disk.sample(n);
        let m = disk.params.mass / n as f64;
        for (r0, r1) in [(1.0, 3.0), (4.0, 6.0), (8.0, 11.0)] {
            let count = parts.iter()
                .filter(|x| (r0..r1)
                .contains(&x[0].hypot(x[1])))
                .count();
            let mass = count as f64 * m;
            let k = 200;
            
            let expected: f64 = (0..k)
                .map(|j| {
                    let r = r0 + (j as f64 + 0.5) * (r1 - r0) / k as f64;
                    2.0 * PI * r * disk.sigma(r) * (r1 - r0) / k as f64
                })
                .sum();
            assert!((mass / expected - 1.0).abs() < 0.01, "annulus {r0}-{r1}: {mass} vs {expected}");
        }
        // Mean |z| in an annulus at R ~ 5 matches the density's vertical profile.
        let ring: Vec<f64> = parts.iter().filter(|x| (4.8..5.2).contains(&x[0].hypot(x[1]))).map(|x| x[2].abs()).collect();
        let mean_z = ring.iter().sum::<f64>() / ring.len() as f64;
        let expected = disk.mean_abs_z(5.0);
        println!("mean |z| at R = 5: sampled {mean_z:.4} kpc, expected {expected:.4} kpc");
        assert!((mean_z / expected - 1.0).abs() < 0.03);
    }

    fn gas_rich() -> GasDiskParams {
        GasDiskParams { mass: 2.0, scale_length: 4.0, r_max: 20.0, sound_speed: 12.0 }
    }

    #[test]
    fn self_gravity_thins_the_disk_and_satisfies_vertical_balance() {
        let mw = milky_way();
        let plain = GasDisk::new(gas_rich(), &mw);
        let opts = DiskOptions { self_gravity: true, ..DiskOptions::default() };
        let sg = GasDisk::with_options(gas_rich(), &mw, opts);
        for r in [2.0, 5.0, 10.0] {
            let (z0, z1) = (plain.mean_abs_z(r), sg.mean_abs_z(r));
            println!("R = {r}: mean |z| {z0:.4} -> {z1:.4} kpc with self-gravity");
            assert!(z1 < z0);
            // Column density is unchanged.
            let top = sg.z_top(r);
            let nz = 20_000;
            let dz = 2.0 * top / nz as f64;
            let col: f64 = (0..nz).map(|k| sg.density(&[r, 0.0, -top + (k as f64 + 0.5) * dz]) * dz).sum();
            assert!((col / sg.sigma(r) - 1.0).abs() < 3e-3, "column at R = {r}");
            // c^2 dln(rho)/dz = -(dPhi_ext/dz + 2 pi G Sigma(<z)) near one scale height.
            let z = z1;
            let h = 1e-3 * z;
            let c2 = 144.0;
            let lhs = c2 * (sg.density(&[r, 0.0, z + h]).ln() - sg.density(&[r, 0.0, z - h]).ln()) / (2.0 * h);
            let n_in = 4000;
            let sigma_in: f64 = (0..n_in).map(|k| 2.0 * sg.density(&[r, 0.0, (k as f64 + 0.5) * z / n_in as f64]) * z / n_in as f64).sum();
            let gz = -mw.acceleration(&[r, 0.0, z], 0.0)[2];
            let rhs = -(gz + 2.0 * PI * G * sigma_in);
            assert!((lhs / rhs - 1.0).abs() < 0.02, "R = {r}: {lhs} vs {rhs}");
        }
    }

    #[test]
    fn noise_displaces_by_the_requested_fraction_of_the_spacing_reproducibly() {
        let mw = milky_way();
        let quiet = GasDisk::new(params(), &mw).sample(5000);
        let opts = DiskOptions { noise_fraction: 0.2, noise_seed: 3, ..DiskOptions::default() };
        let noisy = GasDisk::with_options(params(), &mw, opts);
        let a = noisy.sample(5000);
        let b = noisy.sample(5000);
        assert_eq!(a, b, "same seed, same points");
        let m = params().mass / 5000.0;
        let ratios: Vec<f64> = quiet
            .iter()
            .zip(&a)
            .map(|(x0, x1)| {
                let d = ((x1[0] - x0[0]).powi(2) + (x1[1] - x0[1]).powi(2) + (x1[2] - x0[2]).powi(2)).sqrt();
                d / (m / noisy.density(x0)).cbrt()
            })
            .collect();
        // |3D Gaussian| with sigma = 0.2 per axis has rms 0.2 sqrt(3).
        let rms = (ratios.iter().map(|r| r * r).sum::<f64>() / ratios.len() as f64).sqrt();
        println!("rms displacement / spacing = {rms:.3}");
        assert!((rms / (0.2 * 3f64.sqrt()) - 1.0).abs() < 0.1);
    }
}
