//! Chapter 4, Lesson 4: Low-pass filter driven by white noise
//!
//! Simulates a first-order low-pass filter (time constant TAU) driven by white noise
//! and compares with theoretical standard deviation bounds.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;
use rand_distr::{Distribution, Normal};

/// Simulation results
pub struct Results {
    pub time: Vec<f64>,
    pub y: Vec<f64>,
    pub sig_plus: Vec<f64>,
    pub sig_minus: Vec<f64>,
}

/// Run the C4L4 simulation
pub fn run() -> Results {
    let tau: f64 = 0.2;
    let phi: f64 = 1.0;
    let h: f64 = 0.01;
    let sig = (phi / h).sqrt();

    let normal = Normal::new(0.0, 1.0).unwrap();
    let mut rng = rand::thread_rng();

    let mut y = 0.0;
    let mut t = 0.0;

    let mut array_t = Vec::new();
    let mut array_y = Vec::new();
    let mut array_sig_plus = Vec::new();
    let mut array_sig_minus = Vec::new();

    while t <= 5.0 {
        let x = sig * normal.sample(&mut rng);
        let y_old = y;

        // RK2 integration (second-order Runge-Kutta)
        let yd = (x - y) / tau;
        y += h * yd;
        t += h;
        let yd = (x - y) / tau;
        y = (y_old + y) / 2.0 + 0.5 * h * yd;

        // Theoretical standard deviation bounds
        let sig_plus = (phi * (1.0 - (-2.0_f64 * t / tau).exp()) / (2.0 * tau)).sqrt();
        let sig_minus = -sig_plus;

        array_t.push(t);
        array_y.push(y);
        array_sig_plus.push(sig_plus);
        array_sig_minus.push(sig_minus);
    }

    Results {
        time: array_t,
        y: array_y,
        sig_plus: array_sig_plus,
        sig_minus: array_sig_minus,
    }
}

/// Run and save results to file
pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    // Save data file
    let data_file = format!("{}/c4l4_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.y.clone(),
        results.sig_plus.clone(),
        results.sig_minus.clone(),
    ])?;

    // Create plot
    let plot_file = format!("{}/c4l4_plot.png", output_dir);
    let config = PlotConfig::new("Simulation of low-pass filter driven by white noise")
        .with_labels("Time (S)", "Y");

    let series = vec![
        Series::new(results.time.clone(), results.y.clone())
            .with_label("Y")
            .with_color(plotters::prelude::BLUE),
        Series::new(results.time.clone(), results.sig_plus.clone())
            .with_label("+σ")
            .with_color(plotters::prelude::RED),
        Series::new(results.time.clone(), results.sig_minus.clone())
            .with_label("-σ")
            .with_color(plotters::prelude::RED),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C4L4: Simulation finished");
    println!("  Data saved to: {}", data_file);
    println!("  Plot saved to: {}", plot_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c4l4_runs() {
        let results = run();
        assert!(results.time.len() > 0);
        assert_eq!(results.time.len(), results.y.len());
        assert_eq!(results.time.len(), results.sig_plus.len());
    }
}
