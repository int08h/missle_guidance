//! Chapter 1, Lesson 2: Simple Harmonic Oscillator (Second-Order RK2)
//!
//! Demonstrates second-order Runge-Kutta integration for a simple harmonic oscillator.
//! Shows improved accuracy over first-order Euler.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

/// Simulation results
pub struct Results {
    pub time: Vec<f64>,
    pub x: Vec<f64>,
    pub x_theory: Vec<f64>,
}

/// Run the C1L2 simulation
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

        // Save old values
        let xold = x;
        let xdold = xd;

        // First derivative evaluation
        let xdd = -w * w * x;

        // Euler step
        x += h * xd;
        xd += h * xdd;
        t += h;

        // Second derivative evaluation at new position
        let xdd = -w * w * x;

        // RK2 averaging correction
        x = 0.5 * (xold + x + h * xd);
        xd = 0.5 * (xdold + xd + h * xdd);

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
    let data_file = format!("{}/c1l2_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.x.clone(),
        results.x_theory.clone(),
    ])?;

    // Create plot
    let plot_file = format!("{}/c1l2_plot.png", output_dir);
    let config = PlotConfig::new("Harmonic Oscillator - Second Order RK2")
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

    println!("C1L2: Simulation finished");
    println!("  Data saved to: {}", data_file);
    println!("  Plot saved to: {}", plot_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c1l2_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }

    #[test]
    fn test_c1l2_more_accurate_than_l1() {
        let results = run();
        // RK2 should be more accurate than first-order Euler
        // Check that numerical matches theory closely
        let max_error: f64 = results.x.iter()
            .zip(results.x_theory.iter())
            .map(|(x, xt)| (x - xt).abs())
            .fold(0.0, f64::max);

        assert!(max_error < 0.01, "RK2 should be very accurate, max error: {}", max_error);
    }
}
