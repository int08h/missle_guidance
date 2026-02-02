//! Chapter 3, Lesson 1: PN Miss Distance Analysis
//!
//! Analyzes miss distance contributions from target maneuver and heading error
//! using adjoint method.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

/// Simulation results
pub struct Results {
    pub time: Vec<f64>,
    pub xm_nt: Vec<f64>,   // Miss due to target maneuver
    pub xm_he: Vec<f64>,   // Miss due to heading error
}

/// Run the C3L1 simulation
pub fn run() -> Results {
    let xnt = 96.6;      // Target acceleration (ft/s^2) - 3G
    let xnp = 4.0;       // Navigation ratio
    let tau = 1.0;       // Time constant
    let tf = 10.0;       // Flight time
    let vm = 3000.0;     // Missile velocity
    let he_deg = -20.0;  // Heading error (degrees)
    let h = 0.01;        // Time step

    let t = 0.0;
    let mut s = 0.0;
    let mut tp = t + 0.00001;

    // Adjoint state variables
    let mut x1 = 0.0;
    let mut x2 = 0.0;
    let mut x3 = 1.0;
    let mut x4 = 0.0;

    let he = he_deg / 57.3;

    let mut array_tp = Vec::new();
    let mut array_xm_nt = Vec::new();
    let mut array_xm_he = Vec::new();

    while tp <= tf - 1e-5 {
        // Save old values
        let x1_old = x1;
        let x2_old = x2;
        let x3_old = x3;
        let x4_old = x4;

        // First derivative evaluation
        let x1d = x2;
        let x2d = x3;
        let y1 = (x4 - x2) / tau;
        let tgo = tp + 0.00001;
        let x3d = xnp * y1 / tgo;
        let x4d = -y1;

        // Euler step
        x1 += h * x1d;
        x2 += h * x2d;
        x3 += h * x3d;
        x4 += h * x4d;
        tp += h;

        // Second derivative evaluation
        let x1d = x2;
        let x2d = x3;
        let y1 = (x4 - x2) / tau;
        let tgo = tp + 0.00001;
        let x3d = xnp * y1 / tgo;
        let x4d = -y1;

        // RK2 averaging
        x1 = (x1_old + x1) / 2.0 + 0.5 * h * x1d;
        x2 = (x2_old + x2) / 2.0 + 0.5 * h * x2d;
        x3 = (x3_old + x3) / 2.0 + 0.5 * h * x3d;
        x4 = (x4_old + x4) / 2.0 + 0.5 * h * x4d;

        s += h;

        // Store data at sampling interval
        if s >= 0.0999 {
            s = 0.0;
            array_tp.push(tp);
            array_xm_nt.push(xnt * x1);
            array_xm_he.push(-vm * he * x2);
        }
    }

    Results {
        time: array_tp,
        xm_nt: array_xm_nt,
        xm_he: array_xm_he,
    }
}

/// Run and save results to file
pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    // Save data file
    let data_file = format!("{}/c3l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.xm_nt.clone(),
        results.xm_he.clone(),
    ])?;

    // Create target maneuver miss plot
    let plot_file1 = format!("{}/c3l1_target_miss.png", output_dir);
    let config = PlotConfig::new("Target Maneuver Miss")
        .with_labels("Flight Time (Sec)", "Target Maneuver Miss (Ft)");

    let series = vec![
        Series::new(results.time.clone(), results.xm_nt.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file1, &config, &series).ok();

    // Create heading error miss plot
    let plot_file2 = format!("{}/c3l1_heading_miss.png", output_dir);
    let config = PlotConfig::new("Heading Error Miss")
        .with_labels("Flight Time (Sec)", "Heading Error Miss (Ft)");

    let series = vec![
        Series::new(results.time.clone(), results.xm_he.clone())
            .with_color(plotters::prelude::RED),
    ];

    line_plot(&plot_file2, &config, &series).ok();

    println!("C3L1: Simulation finished");
    println!("  Data saved to: {}", data_file);
    println!("  Plot saved to: {}", plot_file1);
    println!("  Plot saved to: {}", plot_file2);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c3l1_runs() {
        let results = run();
        assert!(!results.time.is_empty());
        assert_eq!(results.time.len(), results.xm_nt.len());
        assert_eq!(results.time.len(), results.xm_he.len());
    }
}
