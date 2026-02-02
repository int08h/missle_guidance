//! Chapter 34, Lesson 2: Adjoint Model Using Shaping Filter Approach
//!
//! Computes RMS miss using adjoint method with shaping filter
//! for target maneuver modeling.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub rms: Vec<f64>,
}

/// Run the C34L2 simulation
pub fn run() -> Results {
    let _xnt: f64 = 96.6;
    let xnp: f64 = 3.0;
    let tau: f64 = 1.0;
    let tf: f64 = 10.0;
    let h: f64 = 0.01;
    let xnu: f64 = 0.5;
    let beta: f64 = 96.6;

    let t: f64 = 0.0;
    let mut s: f64 = 0.0;
    let mut tp = t + 0.00001;

    // State variables
    let mut x1: f64 = 0.0;
    let mut x2: f64 = 0.0;
    let mut x3: f64 = 1.0;
    let mut x4: f64 = 0.0;
    let mut x5: f64 = 0.0;
    let mut x6: f64 = 0.0;
    let mut x7: f64 = 0.0;

    let mut array_tp = Vec::new();
    let mut array_rms = Vec::new();

    while tp <= tf - 1e-5 {
        s += h;
        let x1old = x1;
        let x2old = x2;
        let x3old = x3;
        let x4old = x4;
        let x5old = x5;
        let x6old = x6;
        let x7old = x7;

        // First derivative evaluation
        let x1d = x2;
        let x2d = x3;
        let y1 = (x4 - x2) / tau;
        let tgo = tp + 0.00001;
        let x3d = xnp * y1 / tgo;
        let x4d = -y1;
        let x5d = x1 * x1;
        let x6d = x2 - 2.0 * xnu * x6;
        let x7d = (2.0 * xnu * x6).powi(2);

        // Euler step
        x1 += h * x1d;
        x2 += h * x2d;
        x3 += h * x3d;
        x4 += h * x4d;
        x5 += h * x5d;
        x6 += h * x6d;
        x7 += h * x7d;
        tp += h;

        // Second derivative for RK2
        let x1d = x2;
        let x2d = x3;
        let y1 = (x4 - x2) / tau;
        let tgo = tp + 0.00001;
        let x3d = xnp * y1 / tgo;
        let x4d = -y1;
        let x5d = x1 * x1;
        let x6d = x2 - 2.0 * xnu * x6;
        let x7d = (2.0 * xnu * x6).powi(2);

        // RK2 averaging
        x1 = (x1old + x1) / 2.0 + 0.5 * h * x1d;
        x2 = (x2old + x2) / 2.0 + 0.5 * h * x2d;
        x3 = (x3old + x3) / 2.0 + 0.5 * h * x3d;
        x4 = (x4old + x4) / 2.0 + 0.5 * h * x4d;
        x5 = (x5old + x5) / 2.0 + 0.5 * h * x5d;
        x6 = (x6old + x6) / 2.0 + 0.5 * h * x6d;
        x7 = (x7old + x7) / 2.0 + 0.5 * h * x7d;

        s += h;
        if s >= 0.000999 {
            s = 0.0;
            let xmbt = beta * (x7 / xnu).sqrt();
            let xic = beta * x6;
            let rms = (xmbt.powi(2) + xic.powi(2)).sqrt();

            array_tp.push(tp);
            array_rms.push(rms);
        }
    }

    Results {
        time: array_tp,
        rms: array_rms,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c34l2_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.rms.clone(),
    ])?;

    let plot_file = format!("{}/c34l2_rms.png", output_dir);
    let config = PlotConfig::new("Adjoint Model Using Shaping Filter Approach")
        .with_labels("Flight Time (S)", "RMS Miss (Ft)");

    let series = vec![
        Series::new(results.time.clone(), results.rms.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C34L2: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c34l2_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
