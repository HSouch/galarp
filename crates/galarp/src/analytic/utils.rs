//! General GalaRP based utils, mostly vector operations.
//! 


#[inline]
pub fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
pub fn cross(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}


#[inline]
pub fn norm(a: &[f64; 3]) -> f64 {
    dot(&a, &a).sqrt()
}


pub fn normalized(a: [f64; 3]) -> [f64; 3] {
    let n = norm(&a);
    [a[0] / n, a[1] / n, a[2] / n]
}

/// Unit vectors spanning the plane perpendicular to `los`: x = north x los, y = los x x.
/// If `north` is (anti)parallel to `los`, the coordinate axis least aligned with `los`
/// stands in for it (GalaRP switches to +y in that case).
pub fn basis(los: [f64; 3], north: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let mut x = cross(&north, &los);
    if dot(&x, &x) < 1e-24 * dot(&north, &north).max(f64::MIN_POSITIVE) {
        let k = (0..3).min_by(|&a, &b| los[a].abs().total_cmp(&los[b].abs())).unwrap();
        let mut axis = [0.0; 3];
        axis[k] = 1.0;
        x = cross(&axis, &los);
    }
    let x = normalized(x);
    (x, normalized(cross(&los, &x)))
}