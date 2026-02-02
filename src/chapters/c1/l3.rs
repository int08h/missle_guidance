//! Chapter 1, Lesson 3: Digital Filter Response
//!
//! Demonstrates a first-order digital filter step response and compares
//! with analytical solution.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

/// Simulation results
pub struct Results {
    pub time: Vec<f64>,
    pub y: Vec<f64>,
    pub y_theory: Vec<f64>,
}

/// Run the C1L3 simulation
pub fn run() -> Results {
    let g: f64 = 0.5;   // Filter gain
    let x: f64 = 1.0;   // Step input
    let ts: f64 = 0.1;  // Sample time

    let mut y: f64 = 0.0;

    let mut array_t = Vec::new();
    let mut array_y = Vec::new();
    let mut array_y_theory = Vec::new();

    for n in 1..=20 {
        // Digital filter: y[n] = y[n-1] + G*(x - y[n-1])
        y = y + g * (x - y);
        let t = n as f64 * ts;

        // Theoretical response: y = 1 - (1-G)^n
        let y_theory = 1.0 - (1.0 - g).powi(n);

        array_t.push(t);
        array_y.push(y);
        array_y_theory.push(y_theory);
    }

    Results {
        time: array_t,
        y: array_y,
        y_theory: array_y_theory,
    }
}

/// Run and save results to file
pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    // Save data file
    let data_file = format!("{}/c1l3_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.y.clone(),
        results.y_theory.clone(),
    ])?;

    // Create plot
    let plot_file = format!("{}/c1l3_plot.png", output_dir);
    let config = PlotConfig::new("Digital Filter Step Response")
        .with_labels("T (S)", "Y");

    let series = vec![
        Series::new(results.time.clone(), results.y.clone())
            .with_label("Numerical")
            .with_color(plotters::prelude::BLUE),
        Series::new(results.time.clone(), results.y_theory.clone())
            .with_label("Theory")
            .with_color(plotters::prelude::RED),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C1L3: Simulation finished");
    println!("  Data saved to: {}", data_file);
    println!("  Plot saved to: {}", plot_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c1l3_runs() {
        let results = run();
        assert_eq!(results.time.len(), 20);
    }

    #[test]
    fn test_c1l3_converges_to_one() {
        let results = run();
        let final_y = *results.y.last().unwrap();
        assert!((final_y - 1.0).abs() < 0.001, "Filter should converge to 1.0");
    }

    #[test]
    fn test_c1l3_matches_theory() {
        let results = run();
        for i in 0..results.y.len() {
            let diff = (results.y[i] - results.y_theory[i]).abs();
            assert!(diff < 1e-10, "Numerical should match theory exactly");
        }
    }
}
