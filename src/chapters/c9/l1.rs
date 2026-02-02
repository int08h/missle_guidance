//! Chapter 9, Lesson 1: Kalman Filter Gains
//!
//! Computes position-rate-acceleration Kalman filter gains
//! for a three-state tracking filter.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub k1: Vec<f64>,
    pub k2: Vec<f64>,
    pub k3: Vec<f64>,
}

/// Run the C9L1 simulation
pub fn run() -> Results {
    let vc: f64 = 4000.0;
    let xnt: f64 = 96.6;
    let vm: f64 = 3000.0;
    let hedeg: f64 = 20.0;
    let sigrin: f64 = 0.001;
    let ts: f64 = 0.1;
    let tf: f64 = 10.0;

    let ts2 = ts * ts;
    let ts3 = ts2 * ts;
    let ts4 = ts3 * ts;
    let ts5 = ts4 * ts;

    let phin = xnt * xnt / tf;
    let rtm = vc * tf;
    let sigpos = rtm * sigrin;
    let sign2 = sigpos * sigpos;

    // Covariance matrix elements (symmetric 3x3)
    let mut p11 = sign2;
    let mut p12: f64 = 0.0;
    let mut p13: f64 = 0.0;
    let mut p22 = (vm * hedeg / 57.3).powi(2);
    let mut p23: f64 = 0.0;
    let mut p33 = xnt * xnt;

    let mut t: f64 = 0.0;

    let mut array_t = Vec::new();
    let mut array_k1 = Vec::new();
    let mut array_k2 = Vec::new();
    let mut array_k3 = Vec::new();

    while t <= tf - 1e-5 {
        let tgo = tf - t + 0.000001;
        let rtm = vc * tgo;
        let sigpos = rtm * sigrin;
        let sign2 = sigpos * sigpos;

        // Propagate covariance (M = Phi * P * Phi' + Q)
        let m11 = p11 + ts * p12 + 0.5 * ts2 * p13
            + ts * (p12 + ts * p22 + 0.5 * ts2 * p23)
            + 0.5 * ts2 * (p13 + ts * p23 + 0.5 * ts2 * p33)
            + ts5 * phin / 20.0;
        let m12 = p12 + ts * p22 + 0.5 * ts2 * p23
            + ts * (p13 + ts * p23 + 0.5 * ts2 * p33)
            + ts4 * phin / 8.0;
        let m13 = p13 + ts * p23 + 0.5 * ts2 * p33 + phin * ts3 / 6.0;
        let m22 = p22 + ts * p23 + ts * (p23 + ts * p33) + phin * ts3 / 3.0;
        let m23 = p23 + ts * p33 + 0.5 * ts2 * phin;
        let m33 = p33 + phin * ts;

        // Kalman gains
        let k1 = m11 / (m11 + sign2);
        let k2 = m12 / (m11 + sign2);
        let k3 = m13 / (m11 + sign2);

        // Update covariance
        p11 = (1.0 - k1) * m11;
        p12 = (1.0 - k1) * m12;
        p13 = (1.0 - k1) * m13;
        p22 = -k2 * m12 + m22;
        p23 = -k2 * m13 + m23;
        p33 = -k3 * m13 + m33;

        array_t.push(t);
        array_k1.push(k1);
        array_k2.push(k2);
        array_k3.push(k3);

        t += ts;
    }

    Results {
        time: array_t,
        k1: array_k1,
        k2: array_k2,
        k3: array_k3,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c9l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.k1.clone(),
        results.k2.clone(),
        results.k3.clone(),
    ])?;

    let plot_file = format!("{}/c9l1_plot.png", output_dir);
    let config = PlotConfig::new("Kalman Filter Gains")
        .with_labels("Flight Time (Sec)", "K1");

    let series = vec![
        Series::new(results.time.clone(), results.k1.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C9L1: Simulation finished");
    println!("  Data saved to: {}", data_file);
    println!("  Plot saved to: {}", plot_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c9l1_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }

    #[test]
    fn test_c9l1_gains_bounded() {
        let results = run();
        // K1 should be bounded between 0 and 1
        for k1 in &results.k1 {
            assert!(*k1 >= 0.0 && *k1 <= 1.0);
        }
    }
}
