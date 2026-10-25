
pub mod composite;

pub mod burkert;
pub mod hernquist;
pub mod miyamoto_nagai;
pub mod nfw;
pub mod point_mass;

pub use composite::Composite;

pub use burkert::Burkert;
pub use hernquist::Hernquist;
pub use miyamoto_nagai::MiyamotoNagai;
pub use nfw::Nfw;
pub use point_mass::PointMass;


/// A gravitational potential that may change in time.
pub trait Potential: Send + Sync {
    /// Phi(x, t).
    fn potential(&self, x: &[f64; 3], t: f64) -> f64;

    /// -grad Phi(x, t).
    fn acceleration(&self, x: &[f64; 3], t: f64) -> [f64; 3];
}


/// A potential parameter: constant, or linearly interpolated in time from a table and held
/// at its end values outside the table.
#[derive(Clone, Debug, PartialEq)]
pub enum Param {
    Constant(f64),
    Table { t: Vec<f64>, value: Vec<f64> },
}

impl Param {
    /// Builds a table parameter. `t` must be strictly increasing and the same length as
    /// `value`.
    pub fn table(t: Vec<f64>, value: Vec<f64>) -> Self {
        assert_eq!(t.len(), value.len(), "parameter table columns differ in length");
        assert!(!t.is_empty(), "empty parameter table");
        assert!(t.windows(2).all(|w| w[0] < w[1]), "parameter table times must increase");
        Param::Table { t, value }
    }

    pub fn at(&self, time: f64) -> f64 {
        match self {
            Param::Constant(v) => *v,
            
            Param::Table { t, value } => {
                let k = t.partition_point(|&tk| tk <= time);
                if k == 0 {
                    value[0]
                } else if k == t.len() {
                    value[t.len() - 1]
                } else {
                    let f = (time - t[k - 1]) / (t[k] - t[k - 1]);
                    value[k - 1] + f * (value[k] - value[k - 1])
                }
            }
        }
    }
}

impl From<f64> for Param {
    fn from(v: f64) -> Self {
        Param::Constant(v)
    }
}


#[inline]
pub(crate) fn norm(x: &[f64; 3]) -> f64 {
    (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt()
}


/// Radial acceleration `-g(r) x / r` for a spherical potential, zero at the center.
#[inline]
pub(crate) fn radial(x: &[f64; 3], r: f64, g_over_r: f64) -> [f64; 3] {
    if r == 0.0 {
        return [0.0; 3];
    }
    [-g_over_r * x[0], -g_over_r * x[1], -g_over_r * x[2]]
}


#[cfg(test)]
pub(crate) mod test_util {
    use super::Potential;

    /// Checks the acceleration against a central finite difference of the potential.
    pub fn check_gradient(p: &dyn Potential, points: &[[f64; 3]], t: f64) {
        for x in points {
            let a = p.acceleration(x, t);
            for k in 0..3 {
                let eps = 1e-5 * (1.0 + x[k].abs());
                let mut xp = *x;
                let mut xm = *x;
                xp[k] += eps;
                xm[k] -= eps;
                let fd = -(p.potential(&xp, t) - p.potential(&xm, t)) / (2.0 * eps);
                let scale = a.iter().map(|v| v.abs()).fold(1e-12, f64::max);
                assert!((fd - a[k]).abs() < 1e-6 * scale, "axis {k} at {x:?}: fd {fd}, analytic {}", a[k]);
            }
        }
    }

    /// Positions spanning inner, scale and outer radii, off the symmetry axes.
    pub fn sample_points(scale: f64) -> Vec<[f64; 3]> {
        [0.05, 0.3, 1.0, 3.0, 20.0]
            .iter()
            .map(|&f| [0.6 * f * scale, -0.5 * f * scale, 0.62 * f * scale])
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_parameters_interpolate_and_clamp() {
        let p = Param::table(vec![0.0, 1.0, 3.0], vec![2.0, 4.0, 0.0]);
        assert_eq!(p.at(-1.0), 2.0);
        assert_eq!(p.at(0.5), 3.0);
        assert_eq!(p.at(2.0), 2.0);
        assert_eq!(p.at(5.0), 0.0);
        assert_eq!(Param::from(7.0).at(123.0), 7.0);
    }
}
