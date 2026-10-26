//! Input orbit tables to define galaxy orbits before the run begins
//! 
//! CSV Format: t,x,y,z,vx,vy,vz
//! All in code-units
//! 


use std::fmt;
use std::io::{self, BufRead, Write};


#[derive(Debug)]
pub enum OrbitError {
    Io(io::Error),
    Parse { line: usize, message: String},
    Invalid(String),
}


impl fmt::Display for OrbitError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            OrbitError::Io(e) => write!(f, "{e}"),
            OrbitError::Parse { line, message } => write!(f, "orbit table line {line}: {message}"),
            OrbitError::Invalid(m) => write!(f, "invalid orbit table: {m}"),
        }
    }
}


impl std::error::Error for OrbitError {}


impl From<io::Error> for OrbitError {
    fn from(e: io::Error) -> Self {
        OrbitError::Io(e)
    }
}

const HEADER: &str = "t,x,y,z,vx,vy,vz";


#[derive(Clone, Debug, PartialEq)]
pub struct OrbitTable {
    pub t: Vec<f64>,
    pub pos: Vec<[f64; 3]>,
    pub vel: Vec<[f64; 3]>,
}

impl OrbitTable {

    pub fn new(t: Vec<f64>, pos: Vec<[f64; 3]>, vel: Vec<[f64; 3]>) -> Result<Self, OrbitError> {
        if t.len() < 2 {
            return Err(OrbitError::Invalid("need at least two rows".into()));
        }
        if pos.len() != t.len() || vel.len() != t.len() {
            return Err(OrbitError::Invalid("columns differ in length".into()));
        }
        if let Some(k) = t.windows(2).position(|w| w[0] >= w[1]) {
            return Err(OrbitError::Invalid(format!("times not increasing at row {}", k + 1)));
        }
        Ok(Self { t, pos, vel })
    }

    pub fn t_start(&self) -> f64 {
        self.t[0]
    }

    pub fn t_end(&self) -> f64 {
        self.t[self.t.len() - 1]
    }


    /// Position and velocity at `time`, or `None` outside the table.
    pub fn at(&self, time: f64) -> Option<([f64; 3], [f64; 3])> {
        if !(self.t_start()..=self.t_end()).contains(&time) {
            return None;
        }
        let k = self.t.partition_point(|&tk| tk <= time).clamp(1, self.t.len() - 1);
        let (t0, t1) = (self.t[k - 1], self.t[k]);
        let dt = t1 - t0;
        let s = (time - t0) / dt;
        let (s2, s3) = (s * s, s * s * s);
        let h = [2.0 * s3 - 3.0 * s2 + 1.0, s3 - 2.0 * s2 + s, -2.0 * s3 + 3.0 * s2, s3 - s2];
        let dh = [6.0 * s2 - 6.0 * s, 3.0 * s2 - 4.0 * s + 1.0, -6.0 * s2 + 6.0 * s, 3.0 * s2 - 2.0 * s];
        let (x0, x1, v0, v1) = (self.pos[k - 1], self.pos[k], self.vel[k - 1], self.vel[k]);
        let mut x = [0.0; 3];
        let mut v = [0.0; 3];
        for d in 0..3 {
            x[d] = h[0] * x0[d] + h[1] * dt * v0[d] + h[2] * x1[d] + h[3] * dt * v1[d];
            v[d] = (dh[0] * x0[d] + dh[1] * dt * v0[d] + dh[2] * x1[d] + dh[3] * dt * v1[d]) / dt;
        }
        Some((x, v))
    }

    pub fn write_csv(&self, mut w: impl Write) -> io::Result<()> {
        writeln!(w, "{HEADER}")?;
        for i in 0..self.t.len() {
            let (x, v) = (self.pos[i], self.vel[i]);
            writeln!(
                w,
                "{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
                self.t[i], x[0], x[1], x[2], v[0], v[1], v[2]
            )?;
        }
        Ok(())
    }

    pub fn read_csv(r: impl BufRead) -> Result<Self, OrbitError> {
        let (mut t, mut pos, mut vel) = (Vec::new(), Vec::new(), Vec::new());
        let mut seen_header = false;
        for (n, line) in r.lines().enumerate() {
            let line = line?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            if !seen_header {
                if line.replace(' ', "") != HEADER {
                    return Err(OrbitError::Parse { line: n + 1, message: format!("expected header `{HEADER}`") });
                }
                seen_header = true;
                continue;
            }
            let vals: Vec<f64> = line
                .split(',')
                .map(|s| s.trim().parse::<f64>())
                .collect::<Result<_, _>>()
                .map_err(|e| OrbitError::Parse { line: n + 1, message: e.to_string() })?;
            if vals.len() != 7 {
                return Err(OrbitError::Parse { line: n + 1, message: format!("expected 7 columns, got {}", vals.len()) });
            }
            t.push(vals[0]);
            pos.push([vals[1], vals[2], vals[3]]);
            vel.push([vals[4], vals[5], vals[6]]);
        }
        Self::new(t, pos, vel)
    }

}


#[cfg(test)]
mod tests {
    use super::*;

    /// x(t) = c0 + c1 t + c2 t^2 + c3 t^3 per axis: Hermite interpolation reproduces cubics.
    fn cubic(t: f64) -> ([f64; 3], [f64; 3]) {
        let c = [[1.0, -2.0, 0.5, 0.1], [0.0, 3.0, -1.0, 0.02], [-4.0, 0.0, 0.0, -0.3]];
        let mut x = [0.0; 3];
        let mut v = [0.0; 3];
        for d in 0..3 {
            x[d] = c[d][0] + c[d][1] * t + c[d][2] * t * t + c[d][3] * t * t * t;
            v[d] = c[d][1] + 2.0 * c[d][2] * t + 3.0 * c[d][3] * t * t;
        }
        (x, v)
    }

    fn cubic_table() -> OrbitTable {
        let t: Vec<f64> = vec![0.0, 0.7, 1.5, 3.0];
        let (pos, vel) = t.iter().map(|&s| cubic(s)).unzip();
        OrbitTable::new(t, pos, vel).unwrap()
    }

    #[test]
    fn hermite_interpolation_is_exact_for_cubics() {
        let table = cubic_table();
        for &s in &[0.0, 0.3, 0.7, 1.1, 2.9, 3.0] {
            let (x, v) = table.at(s).unwrap();
            let (xe, ve) = cubic(s);
            for d in 0..3 {
                assert!((x[d] - xe[d]).abs() < 1e-12, "x at {s}");
                assert!((v[d] - ve[d]).abs() < 1e-12, "v at {s}");
            }
        }
        assert!(table.at(-0.1).is_none() && table.at(3.1).is_none());
    }

    #[test]
    fn csv_round_trips_exactly() {
        let table = cubic_table();
        let mut buf = Vec::new();
        table.write_csv(&mut buf).unwrap();
        let back = OrbitTable::read_csv(buf.as_slice()).unwrap();
        assert_eq!(table, back);
    }

    #[test]
    fn rejects_bad_tables() {
        assert!(OrbitTable::new(vec![0.0, 0.0], vec![[0.0; 3]; 2], vec![[0.0; 3]; 2]).is_err());
        assert!(OrbitTable::read_csv("t,x,y,z,vx,vy,vz\n0,1,2\n".as_bytes()).is_err());
    }
}