//! Core utility functions for missile guidance simulations
//!
//! These functions are ports of the MATLAB utility functions from
//! "Tactical and Strategic Missile Guidance"

// Many utility functions are provided for completeness but may not be used by all simulations.
// Allow dead code to avoid warnings for unused helper functions.
#![allow(dead_code)]

pub mod constants;
pub mod lambert3d;
pub mod lambert2d;
pub mod olambert;
pub mod kepler;
pub mod predict;
pub mod project;
pub mod initial;
pub mod gains;
pub mod rk2;

pub use constants::*;
pub use lambert3d::*;
pub use kepler::*;
pub use predict::*;
