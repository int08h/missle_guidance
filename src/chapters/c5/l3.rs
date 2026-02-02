//! Chapter 5, Lesson 3: Adjoint Method with Integration
//!
//! Uses adjoint method to compute RMS acceleration due to target maneuver.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub xmudnt: Vec<f64>,
}

/// Run the C5L3 simulation
pub fn run() -> Results {
    let xnt: f64 = 96.6;      // Target acceleration
    let xnp: f64 = 3.0;       // Navigation ratio
    let tau: f64 = 1.0;       // Time constant
    let tf: f64 = 10.0;       // Flight time
    let tint: f64 = 0.5;      // Integration interval
    let miss: i32 = 0;        // Miss type flag
    let h: f64 = 0.01;        // Time step

    let t: f64 = 0.0;
    let mut s: f64 = 0.0;
    let mut tp: f64 = t + 0.00001 + tint;

    // State variables
    let mut x1: f64 = 0.0;
    let mut x2: f64 = 0.0;
    let mut x3: f64;
    let mut x4: f64;
    let mut x5: f64 = 0.0;

    if miss == 1 {
        x3 = 1.0;
        x4 = 0.0;
    } else {
        x3 = xnp / (tau * tint);
        x4 = -1.0 / tau;
    }

    let mut array_tp = Vec::new();
    let mut array_xmudnt = Vec::new();

    while tp <= tf - 1e-5 {
        let x1old = x1;
        let x2old = x2;
        let x3old = x3;
        let x4old = x4;
        let x5old = x5;

        // First derivative evaluation
        let x1d = x2;
        let x2d = x3;
        let y1 = (x4 - x2) / tau;
        let tgo = tp + 0.00001;
        let x3d = xnp * y1 / tgo;
        let x4d = -y1;
        let x5d = x1 * x1;

        // Euler step
        x1 += h * x1d;
        x2 += h * x2d;
        x3 += h * x3d;
        x4 += h * x4d;
        x5 += h * x5d;
        tp += h;

        // Second derivative for RK2
        let x1d = x2;
        let x2d = x3;
        let y1 = (x4 - x2) / tau;
        let tgo = tp + 0.00001;
        let x3d = xnp * y1 / tgo;
        let x4d = -y1;
        let x5d = x1 * x1;

        // RK2 averaging
        x1 = (x1old + x1) / 2.0 + 0.5 * h * x1d;
        x2 = (x2old + x2) / 2.0 + 0.5 * h * x2d;
        x3 = (x3old + x3) / 2.0 + 0.5 * h * x3d;
        x4 = (x4old + x4) / 2.0 + 0.5 * h * x4d;
        x5 = (x5old + x5) / 2.0 + 0.5 * h * x5d;

        s += h;

        if s >= 0.099999 {
            s = 0.0;
            let tgo = tp + 0.00001;
            let mut xmudnt = xnt * (x5 / tgo).sqrt();
            if miss == 0 {
                xmudnt /= 32.2;
            }
            array_tp.push(tp);
            array_xmudnt.push(xmudnt);
        }
    }

    Results {
        time: array_tp,
        xmudnt: array_xmudnt,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c5l3_datfil.txt", output_dir);
    save_data(&data_file, &[results.time.clone(), results.xmudnt.clone()])?;

    let plot_file = format!("{}/c5l3_plot.png", output_dir);
    let config = PlotConfig::new("Adjoint Method - RMS Acceleration")
        .with_labels("Flight Time (Sec)", "Acceleration (G)");

    let series = vec![
        Series::new(results.time.clone(), results.xmudnt.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C5L3: Simulation finished");
    println!("  Data saved to: {}", data_file);
    println!("  Plot saved to: {}", plot_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c5l3_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
