//! Chapter 7, Lesson 1: Fading Memory Filter
//!
//! Demonstrates noise transmission characteristics of fading memory filter
//! for line-of-sight rate estimation.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;
use rand::prelude::*;
use rand::SeedableRng;
use rand_distr::StandardNormal;

pub struct Results {
    pub time: Vec<f64>,
    pub y: Vec<f64>,
    pub xnc: Vec<f64>,
    pub xlamd: Vec<f64>,
    pub xlamdh: Vec<f64>,
}

/// Run the C7L1 simulation
pub fn run() -> Results {
    run_with_seed(None)
}

/// Run with optional seed for reproducibility
pub fn run_with_seed(seed: Option<u64>) -> Results {
    let vc: f64 = 4000.0;
    let xnt: f64 = 96.6;
    let yic: f64 = 0.0;
    let vm: f64 = 3000.0;
    let hedeg: f64 = 0.0;
    let beta: f64 = 0.3;
    let xnp: f64 = 3.0;
    let signoise: f64 = 0.001;
    let tf: f64 = 10.0;
    let ts: f64 = 0.1;
    let noise: bool = true;
    let h: f64 = 0.01;

    let mut rng: Box<dyn RngCore> = match seed {
        Some(s) => Box::new(rand::rngs::StdRng::seed_from_u64(s)),
        None => Box::new(rand::rngs::StdRng::from_entropy()),
    };

    let mut y: f64 = yic;
    let mut yd: f64 = -vm * hedeg / 57.3;
    let mut t: f64 = 0.0;
    let mut s: f64 = 0.0;

    let gfilter: f64 = 1.0 - beta * beta;
    let hfilter: f64 = (1.0 - beta) * (1.0 - beta);

    let mut xlamh: f64 = 0.0;
    let mut xlamdh: f64 = 0.0;
    let mut xnc: f64 = 0.0;

    let mut array_t = Vec::new();
    let mut array_y = Vec::new();
    let mut array_xnc = Vec::new();
    let mut array_xlamd = Vec::new();
    let mut array_xlamdh = Vec::new();

    while t <= tf - 1e-5 {
        let yold = y;
        let ydold = yd;

        // First derivative evaluation
        let tgo = tf - t + 1e-5;
        let rtm = vc * tgo;
        let xlam = y / (vc * tgo);
        let xlamd = (rtm * yd + y * vc) / (rtm * rtm);
        let ydd = xnt - xnc;

        // Euler step
        y += h * yd;
        yd += h * ydd;
        t += h;

        // Second derivative for RK2
        let ydd = xnt - xnc;

        // RK2 averaging
        y = 0.5 * (yold + y + h * yd);
        yd = 0.5 * (ydold + yd + h * ydd);

        s += h;

        if s >= ts - 1e-5 {
            s = 0.0;

            let xlamnoise = if noise {
                signoise * rng.sample::<f64, _>(StandardNormal)
            } else {
                0.0
            };

            let res = xlam - (xlamh + ts * xlamdh) + xlamnoise;
            xlamh = gfilter * res + xlamh + ts * xlamdh;
            xlamdh += hfilter * res / ts;
            xnc = xnp * vc * xlamdh;

            array_t.push(t);
            array_y.push(y);
            array_xnc.push(xnc);
            array_xlamd.push(xlamd);
            array_xlamdh.push(xlamdh);
        }
    }

    Results {
        time: array_t,
        y: array_y,
        xnc: array_xnc,
        xlamd: array_xlamd,
        xlamdh: array_xlamdh,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c7l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.y.clone(),
        results.xnc.clone(),
        results.xlamd.clone(),
        results.xlamdh.clone(),
    ])?;

    let plot_file = format!("{}/c7l1_plot.png", output_dir);
    let config = PlotConfig::new("Fading Memory Filter - LOS Rate Estimation")
        .with_labels("Time (S)", "Line of Sight Rate (Rad/S)")
        .with_y_range(-0.01, 0.06);

    let series = vec![
        Series::new(results.time.clone(), results.xlamd.clone())
            .with_color(plotters::prelude::BLUE)
            .with_label("Actual"),
        Series::new(results.time.clone(), results.xlamdh.clone())
            .with_color(plotters::prelude::RED)
            .with_label("Estimated"),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C7L1: Simulation finished");
    println!("  Data saved to: {}", data_file);
    println!("  Plot saved to: {}", plot_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c7l1_runs() {
        let results = run_with_seed(Some(12345));
        assert!(!results.time.is_empty());
    }

    #[test]
    fn test_c7l1_filter_tracks() {
        let results = run_with_seed(Some(12345));
        // Filter should produce reasonable estimates during mid-flight
        // (final sample at t=10 has singularity as tgo->0)
        let n = results.xlamd.len();
        if n > 2 {
            let mid_xlamd = results.xlamd[n / 2];
            let mid_xlamdh = results.xlamdh[n / 2];
            // Both values should be reasonably bounded during flight (< 0.1 rad/s)
            assert!(mid_xlamd.abs() < 0.1);
            assert!(mid_xlamdh.abs() < 0.1);
        }
    }
}
