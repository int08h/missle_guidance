//! Missile Guidance Simulations
//!
//! Rust implementation of simulations from "Tactical and Strategic Missile Guidance"
//! by Paul Zarchan.
//!
//! This library provides:
//! - Core utility functions for orbital mechanics and guidance
//! - Chapter-by-chapter simulation implementations
//! - Plotting utilities for visualization

// Allow approximate constants (1.5708, 3.1416, 6.28, etc.) to match MATLAB's exact values.
// This is intentional for numerical compatibility - changing these to std::f64::consts
// would cause output differences that violate the MATLAB matching requirement.
#![allow(clippy::approx_constant)]

pub mod utils;
pub mod chapters;
pub mod plotting;

pub use utils::*;
pub use plotting::*;

/// Save simulation results to a file in ASCII format (matching MATLAB datfil.txt)
pub fn save_data(filename: &str, data: &[Vec<f64>]) -> std::io::Result<()> {
    use std::io::Write;

    let file = std::fs::File::create(filename)?;
    let mut writer = std::io::BufWriter::new(file);

    if data.is_empty() {
        return Ok(());
    }

    let n_rows = data[0].len();

    for i in 0..n_rows {
        for (j, col) in data.iter().enumerate() {
            if j > 0 {
                write!(writer, " ")?;
            }
            write!(writer, "{:14.6e}", col[i])?;
        }
        writeln!(writer)?;
    }

    Ok(())
}

/// Load simulation results from a file
pub fn load_data(filename: &str) -> std::io::Result<Vec<Vec<f64>>> {
    use std::io::BufRead;

    let file = std::fs::File::open(filename)?;
    let reader = std::io::BufReader::new(file);

    let mut data: Vec<Vec<f64>> = Vec::new();
    let mut n_cols = 0;

    for line in reader.lines() {
        let line = line?;
        let values: Vec<f64> = line
            .split_whitespace()
            .filter_map(|s| s.parse().ok())
            .collect();

        if data.is_empty() {
            n_cols = values.len();
            data.resize(n_cols, Vec::new());
        }

        for (i, &val) in values.iter().enumerate() {
            if i < n_cols {
                data[i].push(val);
            }
        }
    }

    Ok(data)
}
