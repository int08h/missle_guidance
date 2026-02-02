//! Chapter 13, Lesson 1: PN Guidance with Dynamics
//!
//! Proportional navigation with missile dynamics (third-order lag)
//! and acceleration limits.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub y: Vec<f64>,
    pub xncg: Vec<f64>,  // Missile acceleration in G
}

/// Run the C13L1 simulation
pub fn run() -> Results {
    let vc: f64 = 4000.0;
    let xnt: f64 = 32.2;         // 1G target acceleration
    let xnclimg: f64 = 7.0;      // Acceleration limit in G
    let yic: f64 = 0.0;
    let vm: f64 = 3000.0;
    let hedeg: f64 = 0.0;
    let tau: f64 = 0.3;          // Time constant
    let xnp: f64 = 3.0;          // Navigation ratio
    let ta: f64 = 5.0;           // Adjoint time
    let r: f64 = -0.01;
    let tf: f64 = 10.0;
    let h: f64 = 0.01;

    let xnclim = xnclimg * 32.2;

    let mut t: f64 = 0.0;
    let mut s: f64 = 0.0;

    let mut y = yic;
    let mut yd = -vm * hedeg / 57.3;
    let mut xnl: f64 = 0.0;
    let mut elamdh: f64 = 0.0;
    let mut x4: f64 = 0.0;
    let mut x5: f64 = 0.0;
    let mut th: f64 = 0.0;
    let mut thh: f64 = 0.0;

    let mut array_t = Vec::new();
    let mut array_y = Vec::new();
    let mut array_xncg = Vec::new();

    while t <= tf - 1e-5 {
        let yold = y;
        let ydold = yd;
        let xnlold = xnl;
        let elamdhhold = elamdh;
        let x4old = x4;
        let x5old = x5;
        let thold = th;
        let thhhold = thh;

        // First derivative evaluation
        let tgo = tf - t + 0.00001;
        let xlam = y / (vc * tgo);
        let eps = xlam - th - thh + r * thh;
        let dd = 5.0 * eps / tau;
        let elamdhd = 5.0 * (dd - elamdh) / tau;

        let mut xnc = xnp * vc * elamdh;
        if xnc > xnclim {
            xnc = xnclim;
        }
        if xnc < -xnclim {
            xnc = -xnclim;
        }

        let x4d = 5.0 * (xnc - x4) / tau;
        let x5d = 5.0 * (x4 - x5) / tau;
        let xnld = 5.0 * (x5 - xnl) / tau;
        let thd = xnl / vm + ta * xnld / vm;
        let thhd = dd - thd;
        let ydd = xnt - xnl;

        // Euler step
        y += h * yd;
        yd += h * ydd;
        xnl += h * xnld;
        elamdh += h * elamdhd;
        x4 += h * x4d;
        x5 += h * x5d;
        th += h * thd;
        thh += h * thhd;
        t += h;

        // Second derivative for RK2
        let tgo = tf - t + 0.00001;
        let xlam = y / (vc * tgo);
        let eps = xlam - th - thh + r * thh;
        let dd = 5.0 * eps / tau;
        let elamdhd = 5.0 * (dd - elamdh) / tau;

        let mut xnc = xnp * vc * elamdh;
        if xnc > xnclim {
            xnc = xnclim;
        }
        if xnc < -xnclim {
            xnc = -xnclim;
        }

        let x4d = 5.0 * (xnc - x4) / tau;
        let x5d = 5.0 * (x4 - x5) / tau;
        let xnld = 5.0 * (x5 - xnl) / tau;
        let thd = xnl / vm + ta * xnld / vm;
        let thhd = dd - thd;
        let ydd = xnt - xnl;

        // RK2 averaging
        y = 0.5 * (yold + y + h * yd);
        yd = 0.5 * (ydold + yd + h * ydd);
        xnl = 0.5 * (xnlold + xnl + h * xnld);
        elamdh = 0.5 * (elamdhhold + elamdh + h * elamdhd);
        x4 = 0.5 * (x4old + x4 + h * x4d);
        x5 = 0.5 * (x5old + x5 + h * x5d);
        th = 0.5 * (thold + th + h * thd);
        thh = 0.5 * (thhhold + thh + h * thhd);

        s += h;

        if s >= 0.0999 {
            s = 0.0;
            array_t.push(t);
            array_y.push(y);
            array_xncg.push(xnc / 32.2);
        }
    }

    Results {
        time: array_t,
        y: array_y,
        xncg: array_xncg,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c13l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.y.clone(),
        results.xncg.clone(),
    ])?;

    let plot_file = format!("{}/c13l1_plot.png", output_dir);
    let config = PlotConfig::new("PN with Dynamics - Missile Acceleration")
        .with_labels("Time (Sec)", "Missile Acceleration (G)");

    let series = vec![
        Series::new(results.time.clone(), results.xncg.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C13L1: Simulation finished");
    println!("  Data saved to: {}", data_file);
    println!("  Plot saved to: {}", plot_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c13l1_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }

    #[test]
    fn test_c13l1_acceleration_limited() {
        let results = run();
        for xncg in &results.xncg {
            assert!(xncg.abs() <= 7.1);  // Within limit (with small tolerance)
        }
    }
}
