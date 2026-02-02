//! Chapter 6, Lesson 2: Higher-Order System with Integrators
//!
//! 10-state adjoint analysis with additional integrator states for
//! step, ramp, and parabolic maneuver miss.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub xmnt: Vec<f64>,
    pub xmntd: Vec<f64>,
    pub xmntdd: Vec<f64>,
}

/// Run the C6L2 simulation
pub fn run() -> Results {
    let xnt: f64 = 32.2;
    let xnp: f64 = 3.0;
    let tau: f64 = 1.0;
    let tf: f64 = 10.0;
    let vc: f64 = 4000.0;
    let xntd: f64 = 32.2;
    let xntdd: f64 = 32.2;
    let h: f64 = 0.01;

    let mut tp: f64 = 0.00001;
    let mut s: f64 = 0.0;

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
    let mut array_xmntd = Vec::new();
    let mut array_xmntdd = Vec::new();

    while tp <= (tf - 1e-5) {
        let x1_old = x1;
        let x2_old = x2;
        let x3_old = x3;
        let x4_old = x4;
        let x5_old = x5;
        let x6_old = x6;
        let x7_old = x7;
        let x8_old = x8;
        let x9_old = x9;
        let x10_old = x10;

        // First derivative evaluation
        let x1d = x2;
        let x2d = x3;
        let y1 = 5.0 * (5.0 * x5 / tau + x4) / tau;
        let tgo = tp + 0.00001;
        let x3d = y1 / (vc * tgo);
        let x4d = -y1;
        let x5d = -5.0 * x5 / tau + 5.0 * x6 * xnp * vc / tau;
        let x6d = -5.0 * x6 / tau + 5.0 * x7 / tau;
        let x7d = -5.0 * x7 / tau + 5.0 * x8 / tau;
        let x8d = -5.0 * x8 / tau - x2;
        let x9d = x1;
        let x10d = x9;

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

        // Second derivative evaluation
        let x1d = x2;
        let x2d = x3;
        let y1 = 5.0 * (5.0 * x5 / tau + x4) / tau;
        let tgo = tp + 0.00001;
        let x3d = y1 / (vc * tgo);
        let x4d = -y1;
        let x5d = -5.0 * x5 / tau + 5.0 * x6 * xnp * vc / tau;
        let x6d = -5.0 * x6 / tau + 5.0 * x7 / tau;
        let x7d = -5.0 * x7 / tau + 5.0 * x8 / tau;
        let x8d = -5.0 * x8 / tau - x2;
        let x9d = x1;
        let x10d = x9;

        // RK2 averaging
        x1 = 0.5 * (x1_old + x1 + h * x1d);
        x2 = 0.5 * (x2_old + x2 + h * x2d);
        x3 = 0.5 * (x3_old + x3 + h * x3d);
        x4 = 0.5 * (x4_old + x4 + h * x4d);
        x5 = 0.5 * (x5_old + x5 + h * x5d);
        x6 = 0.5 * (x6_old + x6 + h * x6d);
        x7 = 0.5 * (x7_old + x7 + h * x7d);
        x8 = 0.5 * (x8_old + x8 + h * x8d);
        x9 = 0.5 * (x9_old + x9 + h * x9d);
        x10 = 0.5 * (x10_old + x10 + h * x10d);

        s += h;
        if s >= 0.0999 {
            s = 0.0;
            array_tp.push(tp);
            array_xmnt.push(xnt * x1);
            array_xmntd.push(xntd * x9);
            array_xmntdd.push(xntdd * x10);
        }
    }

    Results {
        time: array_tp,
        xmnt: array_xmnt,
        xmntd: array_xmntd,
        xmntdd: array_xmntdd,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    // Note: MATLAB saves xmntd twice (bug in original), we match that
    let data_file = format!("{}/c6l2_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.xmnt.clone(),
        results.xmntd.clone(),
        results.xmntd.clone(), // Match MATLAB bug
    ])?;

    let plot_file = format!("{}/c6l2_step.png", output_dir);
    let config = PlotConfig::new("Missile Miss Due To Step Maneuver")
        .with_labels("Normalized Flight Time (Sec)", "Miss (Ft/G-Sec^2)");
    let series = vec![
        Series::new(results.time.clone(), results.xmnt.clone())
            .with_color(plotters::prelude::BLUE),
    ];
    line_plot(&plot_file, &config, &series).ok();

    let plot_file2 = format!("{}/c6l2_ramp.png", output_dir);
    let config2 = PlotConfig::new("Missile Miss Due To Ramp Maneuver")
        .with_labels("Normalized Flight Time (Sec)", "Miss (Ft/G-Sec^3/S)");
    let series2 = vec![
        Series::new(results.time.clone(), results.xmntd.clone())
            .with_color(plotters::prelude::BLUE),
    ];
    line_plot(&plot_file2, &config2, &series2).ok();

    let plot_file3 = format!("{}/c6l2_parabolic.png", output_dir);
    let config3 = PlotConfig::new("Missile Miss Due To Parabolic Maneuver")
        .with_labels("Normalized Flight Time (Sec)", "Miss (Ft/G-Sec^4/S^2)");
    let series3 = vec![
        Series::new(results.time.clone(), results.xmntdd.clone())
            .with_color(plotters::prelude::BLUE),
    ];
    line_plot(&plot_file3, &config3, &series3).ok();

    println!("C6L2: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c6l2_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
