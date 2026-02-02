//! Chapter 38, Lesson 2: Monte Carlo with Trapezoidal Weave Target
//!
//! Monte Carlo simulation with random target maneuver start time
//! using trapezoidal weave pattern.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;
use rand::SeedableRng;
use rand_distr::{Distribution, Uniform};

pub struct Results {
    pub tf: Vec<f64>,
    pub rms: Vec<f64>,
}

/// Compute trapezoidal weave acceleration
fn trapezoidal_weave(t: f64, tstart: f64, xnt: f64, pz: f64, tr: f64, x: f64) -> f64 {
    if t < tstart {
        return 0.0;
    }

    let t_rel = t - tstart;
    let period_num = (t_rel / pz).floor() as i32;
    let tstar = t_rel - (period_num as f64) * pz;

    if tstar < tr / 2.0 {
        2.0 * xnt * tstar / tr
    } else if tstar < tr / 2.0 + x {
        xnt
    } else if tstar < 3.0 * tr / 2.0 + x {
        -2.0 * xnt * tstar / tr + 2.0 * xnt + 2.0 * xnt * x / tr
    } else if tstar < 3.0 * tr / 2.0 + 2.0 * x {
        -xnt
    } else {
        2.0 * xnt * tstar / tr - 4.0 * xnt - 4.0 * xnt * x / tr
    }
}

/// Run the C38L2 simulation
pub fn run() -> Results {
    run_with_seed(None)
}

/// Run with optional seed for reproducibility
pub fn run_with_seed(seed: Option<u64>) -> Results {
    let num_runs: usize = 100;
    let _tstart: f64 = 0.0;
    let pz: f64 = 3.0;
    let tr: f64 = 1.0;
    let _xl = pz / 2.0;
    let vc: f64 = 3000.0;
    let xnt: f64 = 161.0;
    let tau: f64 = 0.5;
    let xnp: f64 = 3.0;
    let _tf: f64 = 10.0;
    let x = pz / 2.0 - tr;

    let mut rng: rand::rngs::StdRng = match seed {
        Some(s) => rand::rngs::StdRng::seed_from_u64(s),
        None => rand::rngs::StdRng::from_entropy(),
    };
    let uniform = Uniform::new(0.0_f64, 1.0_f64);

    let mut array_tf = Vec::new();
    let mut array_rms = Vec::new();

    let mut tf_val = 0.1;
    while tf_val <= 10.0 + 1e-6 {
        let mut z: Vec<f64> = vec![0.0; num_runs];
        let mut z1: f64 = 0.0;

        for (jj, z_item) in z.iter_mut().enumerate().take(num_runs) {
            let sum: f64 = uniform.sample(&mut rng);
            let tstart = sum * tf_val;

            let mut y: f64 = 0.0;
            let mut yd: f64 = 0.0;
            let mut xlamh: f64 = 0.0;
            let mut t: f64 = 0.0;
            let h: f64 = 0.01;
            let mut _s: f64 = 0.0;

            while t <= tf_val - 1e-5 {
                let yold = y;
                let ydold = yd;
                let xlamhold = xlamh;

                // First derivative evaluation
                let ytdd = trapezoidal_weave(t, tstart, xnt, pz, tr, x);
                let tgo = tf_val - t + 0.00001;
                let xlam = y / (vc * tgo);
                let xlamhd = (xlam - xlamh) / tau;
                let xnc = xnp * vc * xlamhd;
                let ydd = ytdd - xnc;

                // Euler step
                y += h * yd;
                yd += h * ydd;
                xlamh += h * xlamhd;
                t += h;

                // Second derivative for RK2
                let ytdd = trapezoidal_weave(t, tstart, xnt, pz, tr, x);
                let tgo = tf_val - t + 0.00001;
                let xlam = y / (vc * tgo);
                let xlamhd = (xlam - xlamh) / tau;
                let xnc = xnp * vc * xlamhd;
                let ydd = ytdd - xnc;

                // RK2 averaging
                y = 0.5 * (yold + y + h * yd);
                yd = 0.5 * (ydold + yd + h * ydd);
                xlamh = 0.5 * (xlamhold + xlamh + h * xlamhd);

                _s += h;
            }

            *z_item = y;
            z1 += *z_item;
            let _xmean = z1 / (jj + 1) as f64;
        }

        // Calculate RMS
        let mut _z1: f64 = 0.0;
        let mut z2: f64 = 0.0;
        let xmean = z.iter().sum::<f64>() / num_runs as f64;

        for z_val in z.iter().take(num_runs) {
            _z1 += (*z_val - xmean).powi(2);
            z2 += z_val.powi(2);
        }

        let rms = if num_runs > 1 {
            (z2 / (num_runs - 1) as f64).sqrt()
        } else {
            0.0
        };

        array_tf.push(tf_val);
        array_rms.push(rms);

        tf_val += 0.1;
    }

    Results {
        tf: array_tf,
        rms: array_rms,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c38l2_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.tf.clone(),
        results.rms.clone(),
    ])?;

    let plot_file = format!("{}/c38l2_rms.png", output_dir);
    let config = PlotConfig::new("RMS Miss for Various Flight Times")
        .with_labels("Flight Time (S)", "RMS MISS (Ft)");

    let series = vec![
        Series::new(results.tf.clone(), results.rms.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C38L2: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c38l2_runs() {
        let results = run_with_seed(Some(12345));
        assert!(!results.tf.is_empty());
    }
}
