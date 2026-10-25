//! Main GalaRP functionality. Essentiall "restricted ram pressure"


pub mod clouds;
pub mod ram_pressure;
pub mod shadow;
pub mod utils;

pub use clouds::Clouds;
pub use ram_pressure::{DensityAt, RamPressure, Wind};
pub use shadow::DynamicShadow;
