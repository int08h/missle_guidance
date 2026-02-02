//! Chapter 4, Lesson 6: Adjoint Model Using Shaping Filter Approach
//!
//! Uses adjoint equations to determine miss distance due to target maneuver
//! as a function of flight time.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

/// Simulation results
pub struct Results {
    pub tp: Vec<f64>,
    pub xmudnt: Vec<f64>,
}

/// Run the C4L6 simulation
pub fn run() -> Results {
    let xnt: f64 = 96.6;
    let xnp: f64 = 3.0;
    let tau: f64 = 1.0;
    let tf: f64 = 10.0;
    let h: f64 = 0.01;

    let t = 0.0;
    let mut s = 0.0;
    let mut tp = t + 0.00001;

    let mut x1 = 0.0;
    let mut x2 = 0.0;
    let mut x3 = 1.0;
    let mut x4 = 0.0;
    let mut x5 = 0.0;

    let mut array_tp = Vec::new();
    let mut array_xmudnt = Vec::new();

    while tp <= (tf - 1e-5) {
        s += h;

        let x1_old = x1;
        let x2_old = x2;
        let x3_old = x3;
        let x4_old = x4;
        let x5_old = x5;

        // First step of RK2
        let x1d = x2;
        let x2d = x3;
        let y1 = (x4 - x2) / tau;
        let tgo = tp + 0.00001;
        let x3d = xnp * y1 / tgo;
        let x4d = -y1;
        let x5d = x1 * x1;

        x1 += h * x1d;
        x2 += h * x2d;
        x3 += h * x3d;
        x4 += h * x4d;
        x5 += h * x5d;
        tp += h;

        // Second step of RK2
        let x1d = x2;
        let x2d = x3;
        let y1 = (x4 - x2) / tau;
        let tgo = tp + 0.00001;
        let x3d = xnp * y1 / tgo;
        let x4d = -y1;
        let x5d = x1 * x1;

        x1 = (x1_old + x1) / 2.0 + 0.5 * h * x1d;
        x2 = (x2_old + x2) / 2.0 + 0.5 * h * x2d;
        x3 = (x3_old + x3) / 2.0 + 0.5 * h * x3d;
        x4 = (x4_old + x4) / 2.0 + 0.5 * h * x4d;
        x5 = (x5_old + x5) / 2.0 + 0.5 * h * x5d;

        s += h;
        if s >= 0.000999 {
            s = 0.0;
            let xmudnt = xnt * (x5 / tgo).sqrt();
            array_tp.push(tp);
            array_xmudnt.push(xmudnt);
        }
    }

    Results {
        tp: array_tp,
        xmudnt: array_xmudnt,
    }
}

/// Run and save results to file
pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    // Save data file
    let data_file = format!("{}/c4l6_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.tp.clone(),
        results.xmudnt.clone(),
    ])?;

    // Create plot
    let plot_file = format!("{}/c4l6_plot.png", output_dir);
    let config = PlotConfig::new("Adjoint model using shaping filter approach")
        .with_labels("Flight Time (S)", "Miss Dist Standard Deviation (Ft)")
        .with_y_range(0.0, 30.0);

    let series = vec![
        Series::new(results.tp.clone(), results.xmudnt.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C4L6: Simulation finished");
    println!("  Data saved to: {}", data_file);
    println!("  Plot saved to: {}", plot_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c4l6_runs() {
        let results = run();
        assert!(results.tp.len() > 0);
        assert_eq!(results.tp.len(), results.xmudnt.len());
    }
}
