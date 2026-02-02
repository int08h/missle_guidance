//! Chapter 28, Lesson 2: Monte Carlo Miss with Shaping Filter
//!
//! Monte Carlo simulation with shaping filter for random target maneuver.
//! Computes RMS miss distance vs flight time.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;
use rand::SeedableRng;
use rand_distr::{Distribution, StandardNormal};

pub struct Results {
    pub tf: Vec<f64>,
    pub rms: Vec<f64>,
}

/// Run the C28L2 simulation
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
    let h: f64 = 0.01;

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

        let phi = beta * beta / tf_val;
        let sig = (phi / h).sqrt();

        for z_item in z.iter_mut() {
            let mut y: f64 = 0.0;
            let mut yd: f64 = 0.0;
            let mut t: f64 = 0.0;
            let mut xnc: f64;
            let mut xnl: f64 = 0.0;
            let mut xnt: f64 = 0.0;

            while t <= tf_val - 1e-5 {
                // Generate random noise
                let randn: f64 = StandardNormal.sample(&mut rng);
                let x: f64 = sig * randn;

                let yold = y;
                let ydold = yd;
                let xnlold = xnl;
                let xntold = xnt;

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
                let ydd = xnt - xnl;
                let xntd = x;

                // Euler step
                y += h * yd;
                yd += h * ydd;
                xnl += h * xnld;
                xnt += h * xntd;
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
                let ydd = xnt - xnl;
                let xntd = x;

                // RK2 averaging
                y = 0.5 * (yold + y + h * yd);
                yd = 0.5 * (ydold + yd + h * ydd);
                xnl = 0.5 * (xnlold + xnl + h * xnld);
                xnt = 0.5 * (xntold + xnt + h * xntd);
            }

            *z_item = y;
            z1 += *z_item;
        }

        let xmean = z1 / num_runs as f64;

        // Calculate RMS
        let mut _z1: f64 = 0.0;
        let mut z2: f64 = 0.0;
        let mut _sigma: f64 = 0.0;
        let mut rms: f64 = 0.0;

        for (i, z_val) in z.iter().enumerate() {
            _z1 += (*z_val - xmean).powi(2);
            z2 += *z_val * *z_val;
            if i == 0 {
                _sigma = 0.0;
                rms = 0.0;
            } else {
                _sigma = (_z1 / i as f64).sqrt();
                rms = (z2 / i as f64).sqrt();
            }
        }

        array_tf.push(tf_val);
        array_rms.push(rms);

        tf_val += 0.2;
    }

    Results {
        tf: array_tf,
        rms: array_rms,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c28l2_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.tf.clone(),
        results.rms.clone(),
    ])?;

    let plot_file = format!("{}/c28l2_rms.png", output_dir);
    let config = PlotConfig::new("RMS Miss for Various Flight Times")
        .with_labels("Flight Time (S)", "RMS Miss (Ft)");

    let series = vec![
        Series::new(results.tf.clone(), results.rms.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C28L2: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c28l2_runs() {
        let results = run_with_seed(Some(12345));
        assert!(!results.tf.is_empty());
    }
}
