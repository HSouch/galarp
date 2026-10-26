//! Main GalaRP functionality. Essentiall "restricted ram pressure"


pub mod clouds;
pub mod integrator;
pub mod ram_pressure;
pub mod shadow;
pub mod utils;

pub use clouds::Clouds;
pub use integrator::{rk4_step, CloudForces};
pub use ram_pressure::{DensityAt, RamPressure, Wind};
pub use shadow::DynamicShadow;
