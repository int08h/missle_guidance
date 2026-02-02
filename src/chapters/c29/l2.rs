//! Chapter 29, Lesson 2: Target Weave Miss Distance vs Flight Time
//!
//! Simulates miss distance for a sinusoidally weaving target using
//! adjoint method with first-order lag flight control system.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub tp: Vec<f64>,
    pub xmweave: Vec<f64>,
}

/// Run the C29L2 simulation
pub fn run() -> Results {
    let xnt: f64 = 193.2;
    let xnp: f64 = 3.0;
    let tau: f64 = 1.0;
    let tf: f64 = 10.0;
    let w: f64 = 3.0;
    let h: f64 = 0.01;

    let t: f64 = 0.0;
    let mut s: f64 = 0.0;
    let mut tp = t + 0.00001;
    let mut x2: f64 = 0.0;
    let mut x3: f64 = 1.0;
    let mut x4: f64 = 0.0;
    let mut x5: f64 = 0.0;
    let mut x6: f64 = 0.0;

    let mut array_tp = Vec::new();
    let mut array_xmweave = Vec::new();

    while tp <= tf - 1e-5 {
        let x2old = x2;
        let x3old = x3;
        let x4old = x4;
        let x5old = x5;
        let x6old = x6;

        // First derivative evaluation
        let x2d = x3;
        let y1 = (-x2 + x4) / tau;
        let tgo = tp + 0.00001;
        let x3d = y1 * xnp / tgo;
        let x4d = -y1;
        let x5d = x2 - w * w * x6;
        let x6d = x5;

        // Euler step
        x2 += h * x2d;
        x3 += h * x3d;
        x4 += h * x4d;
        x5 += h * x5d;
        x6 += h * x6d;
        tp += h;

        // Second derivative evaluation
        let x2d = x3;
        let y1 = (-x2 + x4) / tau;
        let tgo = tp + 0.00001;
        let x3d = y1 * xnp / tgo;
        let x4d = -y1;
        let x5d = x2 - w * w * x6;
        let x6d = x5;

        // RK2 averaging
        x2 = (x2old + x2) / 2.0 + 0.5 * h * x2d;
        x3 = (x3old + x3) / 2.0 + 0.5 * h * x3d;
        x4 = (x4old + x4) / 2.0 + 0.5 * h * x4d;
        x5 = (x5old + x5) / 2.0 + 0.5 * h * x5d;
        x6 = (x6old + x6) / 2.0 + 0.5 * h * x6d;

        s += h;
        if s >= 0.09999 {
            s = 0.0;
            let xmweave = xnt * w * x6;
            array_tp.push(tp);
            array_xmweave.push(xmweave);
        }
    }

    Results {
        tp: array_tp,
        xmweave: array_xmweave,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c29l2_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.tp.clone(),
        results.xmweave.clone(),
    ])?;

    let plot_file = format!("{}/c29l2_miss.png", output_dir);
    let config = PlotConfig::new("Target Weave Miss Distance")
        .with_labels("Flight Time (S)", "Miss (Ft)");

    let series = vec![
        Series::new(results.tp.clone(), results.xmweave.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C29L2: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c29l2_runs() {
        let results = run();
        assert!(!results.tp.is_empty());
    }
}
