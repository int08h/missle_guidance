//! Chapter 8, Lesson 1: Augmented Proportional Navigation
//!
//! Adjoint analysis for Augmented PN (APN) guidance law.
//! Computes miss distance due to target maneuver and heading error.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub xmnt: Vec<f64>,
    pub xmhe: Vec<f64>,
}

/// Run the C8L1 simulation
pub fn run() -> Results {
    let xnt: f64 = 96.6;
    let xnp: f64 = 4.0;
    let tau: f64 = 1.0;
    let tf: f64 = 10.0;
    let vm: f64 = 3000.0;
    let hedeg: f64 = -20.0;
    let apn: f64 = 1.0;  // APN flag (1 = augmented)
    let h: f64 = 0.01;

    let he = hedeg / 57.3;

    let t: f64 = 0.0;
    let mut s: f64 = 0.0;
    let mut tp: f64 = t + 0.00001;

    // State variables
    let mut x1: f64 = 0.0;
    let mut x2: f64 = 0.0;
    let mut x3: f64 = 1.0;
    let mut x4: f64 = 0.0;

    let mut array_tp = Vec::new();
    let mut array_xmnt = Vec::new();
    let mut array_xmhe = Vec::new();

    while tp <= tf - 1e-5 {
        let x1old = x1;
        let x2old = x2;
        let x3old = x3;
        let x4old = x4;

        // First derivative evaluation
        let x1d = x2 + x4 * xnp * apn / (2.0 * tau);
        let x2d = x3 + xnp * x4 / (tau * tp);
        let x3d = xnp * x4 / (tau * tp * tp);
        let x4d = -x4 / tau - x2;

        // Euler step
        x1 += h * x1d;
        x2 += h * x2d;
        x3 += h * x3d;
        x4 += h * x4d;
        tp += h;

        // Second derivative for RK2
        let x1d = x2 + x4 * xnp * apn / (2.0 * tau);
        let x2d = x3 + xnp * x4 / (tau * tp);
        let x3d = xnp * x4 / (tau * tp * tp);
        let x4d = -x4 / tau - x2;

        // RK2 averaging
        x1 = (x1old + x1) / 2.0 + 0.5 * h * x1d;
        x2 = (x2old + x2) / 2.0 + 0.5 * h * x2d;
        x3 = (x3old + x3) / 2.0 + 0.5 * h * x3d;
        x4 = (x4old + x4) / 2.0 + 0.5 * h * x4d;

        s += h;

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
        time: array_tp,
        xmnt: array_xmnt,
        xmhe: array_xmhe,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c8l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.xmnt.clone(),
        results.xmhe.clone(),
    ])?;

    let plot_file = format!("{}/c8l1_plot.png", output_dir);
    let config = PlotConfig::new("Augmented PN - Target Maneuver Miss")
        .with_labels("Flight Time (Sec)", "Target Maneuver Miss (Ft)");

    let series = vec![
        Series::new(results.time.clone(), results.xmnt.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C8L1: Simulation finished");
    println!("  Data saved to: {}", data_file);
    println!("  Plot saved to: {}", plot_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c8l1_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }

    #[test]
    fn test_c8l1_miss_converges() {
        let results = run();
        // Miss should converge towards zero at end of flight
        let final_miss = results.xmnt.last().unwrap().abs();
        assert!(final_miss < 100.0);  // Less than 100 ft miss at end
    }
}
