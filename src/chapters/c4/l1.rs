//! Chapter 4, Lesson 1: Gaussian Random Numbers
//!
//! Generates Gaussian random numbers using the sum of uniform random variables.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;
use rand::Rng;

/// Simulation results
pub struct Results {
    pub index: Vec<f64>,
    pub x: Vec<f64>,
}

/// Run the C4L1 simulation
pub fn run() -> Results {
    let n = 100;
    let mut rng = rand::thread_rng();

    let mut array_i = Vec::new();
    let mut array_x = Vec::new();

    for i in 1..=n {
        // Generate Gaussian using Central Limit Theorem
        // Sum of 12 uniform [0,1) random numbers - 6 gives N(0,1)
        let mut sum = 0.0;
        for _j in 0..12 {
            sum += rng.gen::<f64>();
        }
        let x = sum - 6.0;

        array_i.push(i as f64);
        array_x.push(x);
    }

    Results {
        index: array_i,
        x: array_x,
    }
}

/// Run and save results to file
pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    // Save data file
    let data_file = format!("{}/c4l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.index.clone(),
        results.x.clone(),
    ])?;

    // Create plot
    let plot_file = format!("{}/c4l1_plot.png", output_dir);
    let config = PlotConfig::new("One hundred random numbers with Gaussian distribution")
        .with_labels("Number", "Value");

    let series = vec![
        Series::new(results.index.clone(), results.x.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C4L1: Simulation finished");
    println!("  Data saved to: {}", data_file);
    println!("  Plot saved to: {}", plot_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c4l1_runs() {
        let results = run();
        assert_eq!(results.index.len(), 100);
        assert_eq!(results.x.len(), 100);
    }

    #[test]
    fn test_c4l1_approximately_normal() {
        let results = run();
        let mean: f64 = results.x.iter().sum::<f64>() / results.x.len() as f64;
        // Mean should be close to 0 (within reason for 100 samples)
        assert!(mean.abs() < 0.5, "Mean should be near 0: {}", mean);
    }
}
