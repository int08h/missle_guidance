//! Chapter 42, Lesson 1: Optimal Guidance Monte Carlo
//!
//! Monte Carlo simulation of optimal guidance with Kalman filtering.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;
use rand::Rng;
use rand::SeedableRng;

pub struct Results {
    pub tf: Vec<f64>,
    pub rms: Vec<f64>,
}

/// Run the C42L1 simulation
pub fn run() -> Results {
    run_with_seed(None)
}

/// Run with optional seed for reproducibility
pub fn run_with_seed(seed: Option<u64>) -> Results {
    let vc: f64 = 5.0 * 3280.0;
    let xntic: f64 = 161.0;
    let _vm: f64 = 3000.0;
    let _hedeg: f64 = 20.0;
    let xnp: f64 = 3.0;
    let tau: f64 = 0.2;
    let num_runs: usize = 50;

    let mut rng: rand::rngs::StdRng = match seed {
        Some(s) => rand::rngs::StdRng::seed_from_u64(s),
        None => rand::rngs::StdRng::from_entropy(),
    };

    let mut array_tf = Vec::new();
    let mut array_rms = Vec::new();

    let mut tf = 0.1;
    while tf <= 10.0 {
        let mut z: Vec<f64> = vec![0.0; num_runs];
        let mut _z1: f64 = 0.0;

        for z_item in z.iter_mut() {
            let sum: f64 = rng.gen();
            let tstart = tf * sum;
            let pz: f64 = rng.gen::<f64>() - 0.5;
            let coef = if pz > 0.0 { 1.0 } else { -1.0 };

            let mut y: f64 = 0.0;
            let mut yd: f64 = 0.0;
            let mut t: f64 = 0.0;
            let h: f64 = 0.001;
            let mut xnc: f64 = 0.0;
            let mut xnl: f64 = 0.0;

            while t <= tf - 1e-5 {
                let yold = y;
                let ydold = yd;
                let xnlold = xnl;

                // First derivative
                let xntc = if t < tstart { 0.0 } else { coef * xntic };
                let tgo = tf - t + 0.00001;
                let _rtm = vc * tgo;
                let xnld = (xnc - xnl) / tau;
                let ydd = xntc - xnl;

                // Euler step
                y += h * yd;
                yd += h * ydd;
                xnl += h * xnld;
                t += h;

                // Guidance (simplified)
                let zem = y + yd * tgo;
                xnc = xnp * zem / (tgo * tgo);

                // RK2 averaging
                let xnld = (xnc - xnl) / tau;
                let xntc = if t < tstart { 0.0 } else { coef * xntic };
                let ydd = xntc - xnl;

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
        let rms = if num_runs > 1 {
            (z2 / (num_runs - 1) as f64).sqrt()
        } else {
            0.0
        };

        array_tf.push(tf);
        array_rms.push(rms);

        tf += 0.1;
    }

    Results {
        tf: array_tf,
        rms: array_rms,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c42l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.tf.clone(),
        results.rms.clone(),
    ])?;

    let plot_file = format!("{}/c42l1_rms.png", output_dir);
    let config = PlotConfig::new("Optimal Guidance Monte Carlo")
        .with_labels("Flight Time (S)", "RMS Miss (Ft)");

    let series = vec![
        Series::new(results.tf.clone(), results.rms.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C42L1: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c42l1_runs() {
        let results = run_with_seed(Some(12345));
        assert!(!results.tf.is_empty());
    }
}
