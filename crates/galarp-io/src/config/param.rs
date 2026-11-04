//! `ParamSpec`: a config parameter that is a number or a table in time.

use std::fmt;

use galarp::potentials::Param;
use galarp::units::MYR;
use serde::Deserialize;

use super::{err, ConfigError};

// ---------------------------------------------------------------------------
// Schema
// ---------------------------------------------------------------------------

/// A parameter written either as a number (`a_kpc = 1.0`) or as a table in time
/// (`a_kpc = { t_myr = [0, 100], value = [1.0, 2.0] }`).
#[derive(Clone, Debug, Deserialize)]
#[serde(untagged)]
pub enum ParamSpec {
    Value(f64),
    Table { t_myr: Vec<f64>, value: Vec<f64> },
}

// ---------------------------------------------------------------------------
// Validation
// ---------------------------------------------------------------------------

impl ParamSpec {
    /// Checks a table's shape (`name` is the key, for the message). A number is always fine.
    pub(super) fn validate(&self, name: &str) -> Result<(), ConfigError> {
        if let ParamSpec::Table { t_myr, value } = self {
            if t_myr.is_empty() || t_myr.len() != value.len() {
                return err(format!("{name}: t_myr and value must be non-empty and the same length"));
            }
            if !t_myr.windows(2).all(|w| w[0] < w[1]) {
                return err(format!("{name}: t_myr must increase"));
            }
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Conversion to code units / galarp types
// ---------------------------------------------------------------------------

impl ParamSpec {
    /// Converts to a code-unit `Param`, multiplying values by `scale`.
    pub fn to_param(&self, scale: f64) -> Param {
        match self {
            ParamSpec::Value(v) => Param::Constant(v * scale),
            ParamSpec::Table { t_myr, value } => {
                Param::table(t_myr.iter().map(|t| t * MYR).collect(), value.iter().map(|v| v * scale).collect())
            }
        }
    }
}

impl fmt::Display for ParamSpec {
    /// The value, or "(table in time)".
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParamSpec::Value(v) => write!(f, "{v}"),
            ParamSpec::Table { .. } => f.write_str("(table in time)"),
        }
    }
}
