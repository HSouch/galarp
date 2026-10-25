//! Core library for GalaRP in Bigwig
//! 
//! Gas clouds are test particles in analytic galaxy and cluster potentials, each galaxy
//! carried along a precomputed orbit pushed by an analytic RP wind with column-density
//! shadowing.
//! 

pub mod cluster;
pub mod orbit;
pub mod potentials;

pub mod external;
pub mod orientation;
pub mod rng;
pub mod units;
