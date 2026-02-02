//! Chapter 28, Lesson 1: Monte Carlo Miss Analysis
//!
//! Computes RMS miss distance for various flight times using Monte Carlo
//! simulation with random target maneuvers.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;
use rand::Rng;
use rand::SeedableRng;

pub struct Results {
    pub tf: Vec<f64>,
    pub rms_miss: Vec<f64>,
}

/// Run the C28L1 simulation
pub fn run() -> Results {
    run_with_seed(None)
}

/// Run with optional seed for reproducibility
pub fn run_with_seed(seed: Option<u64>) -> Results {
    let vc: f64 = 4000.0;
    let beta: f64 = 161.0;
    let _vm: f64 = 3000.0;
    let xnp: f64 = 3.0;
    let tau: f64 = 0.2;
    let num_runs: usize = 1000;
    let xlim: f64 = 999999999.0;

    let mut rng: rand::rngs::StdRng = match seed {
        Some(s) => rand::rngs::StdRng::seed_from_u64(s),
        None => rand::rngs::StdRng::from_entropy(),
    };

    let mut array_tf = Vec::new();
    let mut array_rms = Vec::new();

    let mut tf_val = 0.2;
    while tf_val <= 10.0 {
        let mut z: Vec<f64> = vec![0.0; num_runs];
        let mut z1: f64 = 0.0;

        for z_item in z.iter_mut() {
            let sum: f64 = rng.gen();
            let tstart = tf_val * sum;
            let pz: f64 = rng.gen::<f64>() - 0.5;
            let coef = if pz > 0.0 { 1.0 } else { -1.0 };

            let mut y: f64 = 0.0;
            let mut yd: f64 = 0.0;
            let mut t: f64 = 0.0;
            let h: f64 = 0.01;
            let mut xnc: f64;
            let mut xnl: f64 = 0.0;

            while t <= tf_val - 1e-5 {
                let yold = y;
                let ydold = yd;
                let xnlold = xnl;

                // First derivative evaluation
                let tgo = tf_val - t + 0.00001;
                let rtm = vc * tgo;
                let _xlam = y / rtm;
                let xlamd = (rtm * yd + y * vc) / (rtm * rtm);
                xnc = xnp * vc * xlamd;
                if xnc > xlim {
                    xnc = xlim;
                } else if xnc < -xlim {
                    xnc = -xlim;
                }
                let xnld = (xnc - xnl) / tau;
                let xnt_val = if t < tstart { 0.0 } else { coef * beta };
                let ydd = xnt_val - xnl;

                // Euler step
                y += h * yd;
                yd += h * ydd;
                xnl += h * xnld;
                t += h;

                // Second derivative for RK2
                let tgo = tf_val - t + 0.00001;
                let rtm = vc * tgo;
                let _xlam = y / rtm;
                let xlamd = (rtm * yd + y * vc) / (rtm * rtm);
                xnc = xnp * vc * xlamd;
                if xnc > xlim {
                    xnc = xlim;
                } else if xnc < -xlim {
                    xnc = -xlim;
                }
                let xnld = (xnc - xnl) / tau;
                let xnt_val = if t < tstart { 0.0 } else { coef * beta };
                let ydd = xnt_val - xnl;

                // RK2 averaging
                y = 0.5 * (yold + y + h * yd);
                yd = 0.5 * (ydold + yd + h * ydd);
                xnl = 0.5 * (xnlold + xnl + h * xnld);
            }

            *z_item = y;
            z1 += *z_item;
        }

        let _xmean = z1 / num_runs as f64;

        // Calculate RMS
        let mut z2: f64 = 0.0;
        for z_val in z.iter() {
            z2 += *z_val * *z_val;
        }
        let rms = (z2 / (num_runs - 1) as f64).sqrt();

        array_tf.push(tf_val);
        array_rms.push(rms);

        tf_val += 0.2;
    }

    Results {
        tf: array_tf,
        rms_miss: array_rms,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c28l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.tf.clone(),
        results.rms_miss.clone(),
    ])?;

    let plot_file = format!("{}/c28l1_rms.png", output_dir);
    let config = PlotConfig::new("RMS Miss for Various Flight Times")
        .with_labels("Flight Time (S)", "RMS Miss (Ft)");

    let series = vec![
        Series::new(results.tf.clone(), results.rms_miss.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C28L1: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c28l1_runs() {
        let results = run_with_seed(Some(12345));
        assert!(!results.tf.is_empty());
    }
}
