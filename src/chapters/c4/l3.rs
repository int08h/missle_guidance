//! Chapter 4, Lesson 3: Standard Deviation Calculation
//!
//! Demonstrates how sample standard deviation converges as sample size increases.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;
use rand::Rng;

/// Simulation results
pub struct Results {
    pub index: Vec<f64>,
    pub sigma: Vec<f64>,
}

/// Run the C4L3 simulation
pub fn run() -> Results {
    let mut rng = rand::thread_rng();
    let n = 100;

    // Generate random samples
    let mut z = vec![0.0; n];
    let mut z1 = 0.0;

    for z_item in z.iter_mut() {
        let mut sum = 0.0;
        for _j in 0..12 {
            sum += rng.gen::<f64>();
        }
        *z_item = sum - 6.0;
        z1 += *z_item;
    }

    let xmean = z1 / n as f64;

    // Calculate running standard deviation
    let mut array_i = Vec::new();
    let mut array_sigma = Vec::new();

    z1 = 0.0;
    for (i, z_val) in z.iter().enumerate() {
        z1 += (*z_val - xmean).powi(2);
        let sigma = if i == 0 {
            0.0
        } else {
            (z1 / i as f64).sqrt()
        };

        array_i.push((i + 1) as f64);
        array_sigma.push(sigma);
    }

    Results {
        index: array_i,
        sigma: array_sigma,
    }
}

/// Run and save results to file
pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    // Save data file
    let data_file = format!("{}/c4l3_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.index.clone(),
        results.sigma.clone(),
    ])?;

    // Create plot
    let plot_file = format!("{}/c4l3_plot.png", output_dir);
    let config = PlotConfig::new("Sampled Standard Deviation")
        .with_labels("Number of Samples", "Calculated Standard Deviation");

    let series = vec![
        Series::new(results.index.clone(), results.sigma.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C4L3: Simulation finished");
    println!("  Data saved to: {}", data_file);
    println!("  Plot saved to: {}", plot_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c4l3_runs() {
        let results = run();
        assert_eq!(results.index.len(), 100);
        assert_eq!(results.sigma.len(), 100);
    }

    #[test]
    fn test_c4l3_converges_near_one() {
        let results = run();
        let final_sigma = *results.sigma.last().unwrap();
        // For N(0,1), sigma should be near 1.0
        assert!((final_sigma - 1.0).abs() < 0.5, "Sigma should be near 1.0: {}", final_sigma);
    }
}
