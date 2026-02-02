//! Chapter 6, Lesson 1: Higher-Order System Analysis
//!
//! 10-state adjoint analysis for complex flight control system.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub xmnt: Vec<f64>,
}

/// Run the C6L1 simulation
pub fn run() -> Results {
    let qd: f64 = 0.0;
    let xnt: f64 = 32.2;      // 1G target acceleration
    let xnp: f64 = 4.0;
    let t1: f64 = 0.0667;
    let t2: f64 = 0.133;
    let t3: f64 = 0.2;
    let t4: f64 = 0.267;
    let t5: f64 = 0.333;
    let w: f64 = 10.0;
    let z: f64 = 0.7;
    let tf: f64 = 10.0;
    let vc: f64 = 4000.0;
    let h: f64 = 0.01;

    let t: f64 = 0.0;
    let mut s: f64 = 0.0;
    let mut tp: f64 = t + 0.00001;

    // State variables
    let mut x1: f64 = 0.0;
    let mut x2: f64 = 0.0;
    let mut x3: f64 = 1.0;
    let mut x4: f64 = 0.0;
    let mut x5: f64 = 0.0;
    let mut x6: f64 = 0.0;
    let mut x7: f64 = 0.0;
    let mut x8: f64 = 0.0;
    let mut x9: f64 = 0.0;
    let mut x10: f64 = 0.0;

    let mut array_tp = Vec::new();
    let mut array_xmnt = Vec::new();

    while tp <= tf - 1e-5 {
        s += h;

        let x1old = x1;
        let x2old = x2;
        let x3old = x3;
        let x4old = x4;
        let x5old = x5;
        let x6old = x6;
        let x7old = x7;
        let x8old = x8;
        let x9old = x9;
        let x10old = x10;

        // First derivative evaluation
        let tgo = tp + 0.00001;
        let x1d = x2;
        let x2d = x3;
        let x3d = (x4 + x5 / t2) / (vc * tgo * t1);
        let x4d = -(x4 + x5 / t2) / t1;
        let x5d = -x5 / t2 + xnp * vc * x6 / t3;
        let x6d = -x6 / t3 + qd * w * w * x9 + (1.0 - qd) * x7 / t4;
        let x7d = -x7 / t4 + x8 / t5;
        let x8d = -x8 / t5 - x2;
        let x9d = x10 - 2.0 * z * w * x9;
        let x10d = -w * w * x9 - x2;

        // Euler step
        x1 += h * x1d;
        x2 += h * x2d;
        x3 += h * x3d;
        x4 += h * x4d;
        x5 += h * x5d;
        x6 += h * x6d;
        x7 += h * x7d;
        x8 += h * x8d;
        x9 += h * x9d;
        x10 += h * x10d;
        tp += h;

        // Second derivative for RK2
        let tgo = tp + 0.00001;
        let x1d = x2;
        let x2d = x3;
        let x3d = (x4 + x5 / t2) / (vc * tgo * t1);
        let x4d = -(x4 + x5 / t2) / t1;
        let x5d = -x5 / t2 + xnp * vc * x6 / t3;
        let x6d = -x6 / t3 + qd * w * w * x9 + (1.0 - qd) * x7 / t4;
        let x7d = -x7 / t4 + x8 / t5;
        let x8d = -x8 / t5 - x2;
        let x9d = x10 - 2.0 * z * w * x9;
        let x10d = -w * w * x9 - x2;

        // RK2 averaging
        x1 = (x1old + x1) / 2.0 + 0.5 * h * x1d;
        x2 = (x2old + x2) / 2.0 + 0.5 * h * x2d;
        x3 = (x3old + x3) / 2.0 + 0.5 * h * x3d;
        x4 = (x4old + x4) / 2.0 + 0.5 * h * x4d;
        x5 = (x5old + x5) / 2.0 + 0.5 * h * x5d;
        x6 = (x6old + x6) / 2.0 + 0.5 * h * x6d;
        x7 = (x7old + x7) / 2.0 + 0.5 * h * x7d;
        x8 = (x8old + x8) / 2.0 + 0.5 * h * x8d;
        x9 = (x9old + x9) / 2.0 + 0.5 * h * x9d;
        x10 = (x10old + x10) / 2.0 + 0.5 * h * x10d;

        if s >= 0.0999 {
            s = 0.0;
            let xmnt = xnt * x1;
            array_tp.push(tp);
            array_xmnt.push(xmnt);
        }
    }

    Results {
        time: array_tp,
        xmnt: array_xmnt,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c6l1_datfil.txt", output_dir);
    save_data(&data_file, &[results.time.clone(), results.xmnt.clone()])?;

    let plot_file = format!("{}/c6l1_plot.png", output_dir);
    let config = PlotConfig::new("Higher-Order System - Target Maneuver Miss")
        .with_labels("Flight Time (Sec)", "Target Maneuver Miss (Ft)");

    let series = vec![
        Series::new(results.time.clone(), results.xmnt.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C6L1: Simulation finished");
    println!("  Data saved to: {}", data_file);
    println!("  Plot saved to: {}", plot_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c6l1_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
