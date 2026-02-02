//! Chapter 27, Lesson 7: Adjoint Model for Target Maneuver and Heading Error Miss
//!
//! Uses adjoint method to compute target maneuver miss and heading error miss
//! with autopilot time constant TAU.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub tp: Vec<f64>,
    pub xmnt: Vec<f64>,
    pub xmhe: Vec<f64>,
}

/// Run the C27L7 simulation
pub fn run() -> Results {
    let xnp: f64 = 4.0;
    let xnt: f64 = 96.6;
    let tau: f64 = 1.0;
    let tf: f64 = 10.0;
    let vm: f64 = 3000.0;
    let hedeg: f64 = -20.0;
    let h: f64 = 0.01;

    let he = hedeg / 57.3;

    let mut tp = 0.00001;
    let mut s: f64 = 0.0;

    let mut x1: f64 = 0.0;
    let mut x2: f64 = 0.0;
    let mut x3: f64 = 0.0;
    let mut x4: f64 = 0.0;

    let mut array_tp = Vec::new();
    let mut array_xmnt = Vec::new();
    let mut array_xmhe = Vec::new();

    while tp <= tf - 1e-5 {
        s += h;

        let x1old = x1;
        let x2old = x2;
        let x3old = x3;
        let x4old = x4;

        // Impulse at t=0
        let impulse = if tp < h / 2.0 { 2.0 / h } else { 0.0 };

        // First derivative evaluation (FLAG=0)
        let x1d = x2;
        let x2d = x3;
        let y1 = (x4 - x2) / tau;
        let tgo = tp + 0.00001;
        let x3d = xnp * y1 / tgo + impulse;
        let x4d = -y1;

        // Euler step (FLAG=1)
        x1 += h * x1d;
        x2 += h * x2d;
        x3 += h * x3d;
        x4 += h * x4d;
        tp += h;

        // Second derivative for RK2 (impulse=0 in second step)
        let x1d = x2;
        let x2d = x3;
        let y1 = (x4 - x2) / tau;
        let tgo = tp + 0.00001;
        let x3d = xnp * y1 / tgo;
        let x4d = -y1;

        // RK2 averaging
        x1 = (x1old + x1) / 2.0 + 0.5 * h * x1d;
        x2 = (x2old + x2) / 2.0 + 0.5 * h * x2d;
        x3 = (x3old + x3) / 2.0 + 0.5 * h * x3d;
        x4 = (x4old + x4) / 2.0 + 0.5 * h * x4d;

        if s >= 0.0999 {
            s = 0.0;
            let xmnt = xnt * x1;
            let xmhe = -vm * he * x2;
            array_tp.push(tp);
            array_xmnt.push(xmnt);
            array_xmhe.push(xmhe);
        }
    }

    Results {
        tp: array_tp,
        xmnt: array_xmnt,
        xmhe: array_xmhe,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c27l7_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.tp.clone(),
        results.xmnt.clone(),
        results.xmhe.clone(),
    ])?;

    let plot_file1 = format!("{}/c27l7_miss_target.png", output_dir);
    let config1 = PlotConfig::new("Target Maneuver Miss")
        .with_labels("Flight Time (Sec)", "Target Maneuver Miss (Ft)");

    let series1 = vec![
        Series::new(results.tp.clone(), results.xmnt.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file1, &config1, &series1).ok();

    let plot_file2 = format!("{}/c27l7_miss_heading.png", output_dir);
    let config2 = PlotConfig::new("Heading Error Miss")
        .with_labels("Flight Time (Sec)", "Heading Error Miss (Ft)");

    let series2 = vec![
        Series::new(results.tp.clone(), results.xmhe.clone())
            .with_color(plotters::prelude::RED),
    ];

    line_plot(&plot_file2, &config2, &series2).ok();

    println!("C27L7: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c27l7_runs() {
        let results = run();
        assert!(!results.tp.is_empty());
    }
}
