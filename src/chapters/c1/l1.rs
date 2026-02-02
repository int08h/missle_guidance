//! Chapter 1, Lesson 1: Simple Harmonic Oscillator (First-Order Euler)
//!
//! Demonstrates first-order Euler integration for a simple harmonic oscillator.
//! Compares numerical solution with analytical solution.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

/// Simulation results
pub struct Results {
    pub time: Vec<f64>,
    pub x: Vec<f64>,
    pub x_theory: Vec<f64>,
}

/// Run the C1L1 simulation
pub fn run() -> Results {
    let w: f64 = 2.0;  // Angular frequency
    let h: f64 = 0.01; // Time step

    let mut t: f64 = 0.0;
    let mut s: f64 = 0.0;
    let mut x: f64 = 0.0;
    let mut xd: f64 = w;

    let mut array_t = Vec::new();
    let mut array_x = Vec::new();
    let mut array_x_theory = Vec::new();

    while t <= 10.0 {
        s += h;

        // Simple Euler integration (first-order)
        let xdd = -w * w * x;
        xd += h * xdd;
        x += h * xd;
        t += h;

        // Store data at sampling interval
        if s >= 0.09999 {
            s = 0.0;
            let x_theory = (w * t).sin();
            array_t.push(t);
            array_x.push(x);
            array_x_theory.push(x_theory);
        }
    }

    Results {
        time: array_t,
        x: array_x,
        x_theory: array_x_theory,
    }
}

/// Run and save results to file
pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    // Save data file
    let data_file = format!("{}/c1l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.x.clone(),
        results.x_theory.clone(),
    ])?;

    // Create plot
    let plot_file = format!("{}/c1l1_plot.png", output_dir);
    let config = PlotConfig::new("Harmonic Oscillator - First Order Euler")
        .with_labels("Time (Sec)", "X & XTHEORY")
        .with_x_range(0.0, 10.0)
        .with_y_range(-2.0, 2.0);

    let series = vec![
        Series::new(results.time.clone(), results.x.clone())
            .with_label("Numerical")
            .with_color(plotters::prelude::BLUE),
        Series::new(results.time.clone(), results.x_theory.clone())
            .with_label("Theory")
            .with_color(plotters::prelude::RED),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C1L1: Simulation finished");
    println!("  Data saved to: {}", data_file);
    println!("  Plot saved to: {}", plot_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c1l1_runs() {
        let results = run();
        assert!(!results.time.is_empty());
        assert_eq!(results.time.len(), results.x.len());
        assert_eq!(results.time.len(), results.x_theory.len());
    }

    #[test]
    fn test_c1l1_oscillates() {
        let results = run();
        // Check that both numerical and theory oscillate
        let x_min = results.x.iter().cloned().fold(f64::INFINITY, f64::min);
        let x_max = results.x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(x_min < 0.0);
        assert!(x_max > 0.0);
    }
}
