//! Random numbers
//! Uses a wrapper around rand_xoshiro
//! 


use rand::RngExt;
use rand_distr::Distribution;
use rand_xoshiro::rand_core::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;


#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Rng(Xoshiro256StarStar);


impl Rng {
    pub fn new(seed: u64) -> Self {
        Self(Xoshiro256StarStar::seed_from_u64(seed))
    }

    pub fn uniform(&mut self) -> f64 {
        self.0.random()
    }

    pub fn sample<T, D: Distribution<T>>(&mut self, dist: &D) -> T {
        dist.sample(&mut self.0)
    }
}


#[cfg(test)]
mod tests {
    use std::{assert_eq, assert_ne};

use super::*;

    #[test]
    fn uniform_vals_in_range() {
        let mut a = Rng::new(42);

        let n = 1000;

        for _ in 0..n {
            let x = a.uniform();
            assert!((0.0 .. 1.0).contains(&x));
        }
    }

    #[test]
    fn rng_is_reproducible() {
        let mut a = Rng::new(42);
        let mut b = Rng::new(42);
        let n = 100_000;
        
        for _ in 0..n {
            let x = a.uniform();
            
            assert_eq!(x, b.uniform());
        }
    }

    #[test]
    fn rng_has_right_mean() {
        let mut a = Rng::new(42);

        let n = 100_000;

        let mut sum = 0.0;      

        for _ in 0..n {
            let x = a.uniform();
            sum += x;
        }
        assert!((sum / n as f64 - 0.5).abs() < 0.01);
    }

    #[test]
    fn rng_diff_seeds_diff_outputs() {
        assert_ne!(Rng::new(1).uniform(), Rng::new(2).uniform());
    }
}