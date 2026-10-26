

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Orientation {
    m: [[f64; 3]; 3],          // Rotation matrix
}

impl Default for Orientation {
    fn default() -> Self {
        Self::IDENTITY
    }
}


impl Orientation {
    pub const IDENTITY: Self = Self {
        m: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    };

    pub fn from_normal(normal: [f64; 3]) -> Self {
        let len = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
        assert!(len > 0.0, "disk normal must be nonzero");
        
        let n = [normal[0] / len, normal[1] / len, normal[2] / len];
        let c = n[2]; // cos of the angle between +z and n
        if c < -1.0 + 1e-12 {
            // Upside down: turn by pi about x.
            return Self { m: [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]] };
        }

        let (vx, vy) = (-n[1], n[0]);
        let f = 1.0 / (1.0 + c);
        let m = [
            [1.0 - vy * vy * f, vx * vy * f, vy],
            [vx * vy * f, 1.0 - vx * vx * f, -vx],
            [-vy, vx, 1.0 - (vx * vx + vy * vy) * f],
        ];

        Self { m }

    
    }

    pub fn from_angles(inclination_deg: f64, azimuth_deg: f64) -> Self {
        let (i, a) = (inclination_deg.to_radians(), azimuth_deg.to_radians());
        Self::from_normal([i.sin() * a.cos(), i.sin() * a.sin(), i.cos()])
    }


    /// The disk's spin axis in the cluster frame.
    pub fn normal(&self) -> [f64; 3] {
        self.to_cluster(&[0.0, 0.0, 1.0])
    }

    pub fn is_identity(&self) -> bool {
        *self == Self::IDENTITY
    }


    /// Galaxy-frame vector to the cluster frame (m v).
    #[inline]
    pub fn to_cluster(&self, v: &[f64; 3]) -> [f64; 3] {
        let m = &self.m;
        [
            m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
            m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
            m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
        ]
    }

    /// Cluster-frame vector to the galaxy frame (m^T v).
    #[inline]
    pub fn to_galaxy(&self, v: &[f64; 3]) -> [f64; 3] {
        let m = &self.m;
        [
            m[0][0] * v[0] + m[1][0] * v[1] + m[2][0] * v[2],
            m[0][1] * v[0] + m[1][1] * v[1] + m[2][1] * v[2],
            m[0][2] * v[0] + m[1][2] * v[1] + m[2][2] * v[2],
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: [f64; 3], b: [f64; 3]) -> bool {
        (0..3).all(|k| (a[k] - b[k]).abs() < 1e-12)
    }

    #[test]
    fn rotations_take_z_to_the_normal_and_are_orthonormal() {
        for n in [[1.0, 0.0, 0.0], [0.0, 1.0, 1.0], [-0.3, 0.2, -0.9], [0.0, 0.0, -2.0], [0.0, 0.0, 5.0], [1e-9, 0.0, 1.0]] {
            let o = Orientation::from_normal(n);
            let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            assert!(close(o.normal(), [n[0] / len, n[1] / len, n[2] / len]), "{n:?}");
            for v in [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.3, -0.4, 0.5]] {
                assert!(close(o.to_galaxy(&o.to_cluster(&v)), v));
                let w = o.to_cluster(&v);
                let norm = |a: [f64; 3]| (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
                assert!((norm(w) - norm(v)).abs() < 1e-12);
            }
            // Right-handed: x' cross y' = z'.
            let (x, y, z) = (o.to_cluster(&[1.0, 0.0, 0.0]), o.to_cluster(&[0.0, 1.0, 0.0]), o.normal());
            assert!(close([x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0]], z));
        }
        assert!(Orientation::from_normal([0.0, 0.0, 3.0]).is_identity());
    }

    #[test]
    fn angles_give_the_expected_normal() {
        assert!(close(Orientation::from_angles(90.0, 0.0).normal(), [1.0, 0.0, 0.0]));
        assert!(close(Orientation::from_angles(90.0, 90.0).normal(), [0.0, 1.0, 0.0]));
        let s = 0.5f64.sqrt();
        assert!(close(Orientation::from_angles(45.0, 180.0).normal(), [-s, 0.0, s]));
    }
}
