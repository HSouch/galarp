//! CSV output logs for diagnostic stuff.
//! 


use std::fs::{File, OpenOptions};
use std::io::{self, BufWriter, Write};
use std::path::Path;


pub struct CsvLog {
    out: BufWriter<File>,
    columns: usize,
}


impl CsvLog {
    /// Make a new CSV log file
    pub fn create(path: impl AsRef<Path>, header: &[&str]) -> io::Result<Self> {
        let mut out = BufWriter::new(File::create(path)?);

        writeln!(out, "{}", header.join(","))?;

        out.flush()?;
        Ok(Self { out, columns: header.len() })
    }

    /// Continue an existing log post-restart.
    /// Keeps header and rows before the restart time (t_max)
    /// Creates file if it doesn't existe, fails if existing header has different column names
    pub fn resume(path: impl AsRef<Path>, header: &[&str], t_max: f64) -> io::Result<Self> {

        let path = path.as_ref();
        let Ok(text) = std::fs::read_to_string(path) else {
            return Self::create(path, header);
        };

        let mut lines = text.lines();
        let expected = header.join(",");

        if lines.next() != Some(expected.as_str()) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData, 
                format!("{}: header differs from this run's", path.display())
            ));
        }
        
        let tol = 1e-9 * t_max.abs().max(1.0);
        let kept: Vec<&str> = lines
            .filter(|l| l.split(',').next()
            .and_then(|t| t.parse::<f64>().ok())
            .is_some_and(|t| t <= t_max + tol)).collect();

        let mut out = BufWriter::new(File::create(path)?);

        writeln!(out, "{expected}")?;
        for l in kept {
            writeln!(out, "{l}")?;
        } 

        out.flush()?;
        drop(out);

        let out = BufWriter::new(OpenOptions::new().append(true).open(path)?);
        Ok(Self { out, columns: header.len() })
    }
    
    /// Generate row from a set of values.
    pub fn row(&mut self, values: &[f64]) -> io::Result<()> {
        assert_eq!(values.len(), self.columns, "row length does not match header");

        let line: Vec<String> = values.iter()
            .map(|v| format!("{v:.10e}"))
            .collect();

        writeln!(self.out, "{}", line.join(","))?;
        self.out.flush()
    }
}


