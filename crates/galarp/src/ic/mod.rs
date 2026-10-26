

pub mod disk;

pub use disk::{DiskOptions, GasDisk, GasDiskParams};


pub fn halton(mut index: u64, base: u64) -> f64 {
    let mut f = 1.0;
    let mut r = 0.0;
    let b = base as f64;

    while index > 0 {
        f /= b;
        r += f * (index % base) as f64;
        index /= base;
    }
    
    r

}


#[cfg(test)]

mod tests {
    use super::*;

    #[test]
    fn halton_matches_known_values() {
        assert_eq!(halton(1, 2), 0.5);
        assert_eq!(halton(2, 2), 0.25);
        assert_eq!(halton(3, 2), 0.75);

        assert!((halton(1, 3) - 1.0 / 3.0).abs() < 1e-15);
        assert!((halton(4, 3) - 4.0 / 9.0).abs() < 1e-15);
    }

    #[test]
    fn halton_fills_evenly() {
        let mut bins = [0; 100];

        for i in 0..1000 {
            bins[(halton(i, 10) * 100.0 + 1e-9 ) as usize] += 1;
        }

        assert!(bins.iter().all(|&b| b == 10));
    }
}