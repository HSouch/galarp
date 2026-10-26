//! Sum of potential terms, e.g. disk + bulge + halo.

use super::Potential;

#[derive(Default)]
pub struct Composite {
    pub terms: Vec<Box<dyn Potential>>,
}

impl Composite {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with(mut self, term: impl Potential + 'static) -> Self {
        self.terms.push(Box::new(term));
        self
    }
}

impl Potential for Composite {
    fn potential(&self, x: &[f64; 3], t: f64) -> f64 {
        self.terms.iter().map(|p| p.potential(x, t)).sum()
    }

    fn acceleration(&self, x: &[f64; 3], t: f64) -> [f64; 3] {
        let mut a = [0.0; 3];
        for p in &self.terms {
            let ap = p.acceleration(x, t);
            for k in 0..3 {
                a[k] += ap[k];
            }
        }
        a
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::potentials::test_util::{check_gradient, sample_points};
    use crate::potentials::{Hernquist, MiyamotoNagai, Nfw, Param};

    #[test]
    fn galaxy_model_sums_terms_and_follows_time_dependent_parameters() {
        let galaxy = Composite::new()
            .with(MiyamotoNagai { mass: 5.0.into(), a: 3.0.into(), b: 0.3.into() })
            .with(Hernquist { mass: 1.0.into(), a: 0.5.into() })
            .with(Nfw { m: Param::table(vec![0.0, 1.0], vec![60.0, 30.0]), r_s: 15.0.into() });
        check_gradient(&galaxy, &sample_points(3.0), 0.5);

        // Halving the halo's scale mass weakens the far field.
        let far = [500.0, 0.0, 0.0];
        assert!(galaxy.acceleration(&far, 1.0)[0].abs() < galaxy.acceleration(&far, 0.0)[0].abs());
    }
}
