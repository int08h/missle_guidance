//! Chapter 4, Lesson 2: Gaussian PDF Histogram
//!
//! Creates a histogram of Gaussian random numbers and compares with theoretical PDF.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;
use rand::Rng;
use std::f64::consts::PI;

/// Simulation results
pub struct Results {
    pub ab: Vec<f64>,      // Bin centers
    pub pdf: Vec<f64>,     // Sampled PDF
    pub theory: Vec<f64>,  // Theoretical PDF
}

/// Run the C4L2 simulation
pub fn run() -> Results {
    let xmax = 6.0;
    let xmin = -6.0;
    let range = xmax - xmin;
    let tmp = 1.0 / (2.0 * PI).sqrt();
    let n = 100;
    let bin = 50;

    let mut rng = rand::thread_rng();

    // Generate random samples
    let mut x = vec![0.0; n];
    for x_item in x.iter_mut() {
        let mut sum = 0.0;
        for _j in 0..12 {
            sum += rng.gen::<f64>();
        }
        *x_item = sum - 6.0;
    }

    // Build histogram
    let mut h = vec![0.0; bin];
    for x_val in x.iter() {
        let mut k = (((*x_val - xmin) / range) * bin as f64) as usize;
        if k < 1 { k = 0; }
        if k >= bin { k = bin - 1; }
        h[k] += 1.0;
    }

    // Calculate PDF and theoretical values
    let mut array_ab = Vec::new();
    let mut array_pdf = Vec::new();
    let mut array_th = Vec::new();

    for (k, h_val) in h.iter().enumerate() {
        let pdf = (*h_val / n as f64) * bin as f64 / range;
        let ab = xmin + (k as f64 + 1.0) * range / bin as f64;
        let th = tmp * (-ab * ab / 2.0).exp();

        array_ab.push(ab);
        array_pdf.push(pdf);
        array_th.push(th);
    }

    Results {
        ab: array_ab,
        pdf: array_pdf,
        theory: array_th,
    }
}

/// Run and save results to file
pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    // Save data file
    let data_file = format!("{}/c4l2_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.ab.clone(),
        results.pdf.clone(),
        results.theory.clone(),
    ])?;

    // Create plot
    let plot_file = format!("{}/c4l2_plot.png", output_dir);
    let config = PlotConfig::new("Sample Gaussian distribution")
        .with_labels("X", "Probability Density Function");

    let series = vec![
        Series::new(results.ab.clone(), results.pdf.clone())
            .with_label("Sampled")
            .with_color(plotters::prelude::BLUE),
        Series::new(results.ab.clone(), results.theory.clone())
            .with_label("Theory")
            .with_color(plotters::prelude::RED),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C4L2: Simulation finished");
    println!("  Data saved to: {}", data_file);
    println!("  Plot saved to: {}", plot_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c4l2_runs() {
        let results = run();
        assert_eq!(results.ab.len(), 50);
        assert_eq!(results.pdf.len(), 50);
        assert_eq!(results.theory.len(), 50);
    }
}
