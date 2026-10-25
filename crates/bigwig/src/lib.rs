//! Core library for GalaRP in Bigwig
//! 
//! Gas clouds are test particles in analytic galaxy and cluster potentials, each galaxy
//! carried along a precomputed orbit pushed by an analytic RP wind with column-density
//! shadowing.
//! 

// Larger modules
pub mod cluster;
pub mod galarp;
pub mod ic;
pub mod orbit;
pub mod potentials;

// Single file crates
pub mod external;
pub mod orientation;
pub mod rng;
pub mod units;
