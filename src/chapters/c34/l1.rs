//! Chapter 34, Lesson 1: Shaping Filter Monte Carlo
//!
//! Monte Carlo simulation with random target maneuver timing based on
//! an exponential distribution (shaping filter).

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;
use rand::Rng;
use rand::SeedableRng;
use rand_distr::{Distribution, Normal};

pub struct Results {
    pub tf: Vec<f64>,
    pub rms_miss: Vec<f64>,
}

/// Run the C34L1 simulation
pub fn run() -> Results {
    run_with_seed(None)
}

/// Run with optional seed for reproducibility
pub fn run_with_seed(seed: Option<u64>) -> Results {
    let vc: f64 = 4000.0;
    let _vm: f64 = 3000.0;
    let xnp: f64 = 3.0;
    let tau: f64 = 1.0;
    let num_runs: usize = 1000;
    let beta: f64 = 96.6;
    let xnu: f64 = 0.5;

    let mut rng: rand::rngs::StdRng = match seed {
        Some(s) => rand::rngs::StdRng::seed_from_u64(s),
        None => rand::rngs::StdRng::from_entropy(),
    };
    let normal = Normal::new(0.0, 1.0).unwrap();

    let mut array_tf = Vec::new();
    let mut array_rms = Vec::new();

    let mut tf_val = 0.2;
    while tf_val <= 10.0 {
        let mut z: Vec<f64> = vec![0.0; num_runs];
        let mut _z1: f64 = 0.0;

        for z_item in z.iter_mut() {
            let sig = 1.0 / (2.0 * xnu).sqrt();
            let pz: f64 = rng.gen::<f64>() - 0.5;
            let coef = if pz > 0.0 { 1.0 } else { -1.0 };
            let mut xnt = coef * beta;

            // Initialize random maneuver timing
            let mut qfirst = true;
            let mut delt = 9999.0;
            let mut tnow = 0.0;

            let mut y: f64 = 0.0;
            let mut yd: f64 = 0.0;
            let mut t: f64 = 0.0;
            let h: f64 = 0.01;
            let mut xnc: f64;
            let mut xnl: f64 = 0.0;

            while t <= tf_val - 1e-5 {
                // Shaping filter maneuver timing
                if qfirst {
                    let xnoise1 = sig * normal.sample(&mut rng);
                    let xnoise2 = sig * normal.sample(&mut rng);
                    delt = xnoise1 * xnoise1 + xnoise2 * xnoise2;
                    qfirst = false;
                    tnow = t;
                }
                if t >= delt + tnow {
                    xnt = -xnt;
                    qfirst = true;
                }

                let yold = y;
                let ydold = yd;
                let xnlold = xnl;

                // First derivative evaluation
                let tgo = tf_val - t + 0.00001;
                let rtm = vc * tgo;
                let xlamd = (rtm * yd + y * vc) / (rtm * rtm);
                xnc = xnp * vc * xlamd;
                let xnld = (xnc - xnl) / tau;
                let ydd = xnt - xnl;

                // Euler step
                y += h * yd;
                yd += h * ydd;
                xnl += h * xnld;
                t += h;

                // Second derivative for RK2
                let tgo = tf_val - t + 0.00001;
                let rtm = vc * tgo;
                let xlamd = (rtm * yd + y * vc) / (rtm * rtm);
                xnc = xnp * vc * xlamd;
                let xnld = (xnc - xnl) / tau;
                let ydd = xnt - xnl;

                // RK2 averaging
                y = 0.5 * (yold + y + h * yd);
                yd = 0.5 * (ydold + yd + h * ydd);
                xnl = 0.5 * (xnlold + xnl + h * xnld);
            }

            *z_item = y;
            _z1 += *z_item;
        }

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

    let data_file = format!("{}/c34l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.tf.clone(),
        results.rms_miss.clone(),
    ])?;

    let plot_file = format!("{}/c34l1_rms.png", output_dir);
    let config = PlotConfig::new("Shaping Filter Monte Carlo - RMS Miss")
        .with_labels("Time", "RMS Miss (ft)");

    let series = vec![
        Series::new(results.tf.clone(), results.rms_miss.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C34L1: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c34l1_runs() {
        let results = run_with_seed(Some(12345));
        assert!(!results.tf.is_empty());
    }
}
