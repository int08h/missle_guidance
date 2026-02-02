//! Chapter 7, Lesson 2: Digital Fading Memory Filter Monte Carlo
//!
//! Monte Carlo simulation with fading memory filter for noise analysis.
//!
//! MATLAB BUG: The original MATLAB code calls `gaussc7(SIGNOISE)` which is an
//! undefined function. We fix this by using `SIGNOISE * randn` as noted in CLAUDE.md.
//! Output: 3 columns [TF, SIGMA, XMEAN]

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;
use rand_distr::{Distribution, Normal};

pub struct Results {
    pub tf: Vec<f64>,
    pub sigma: Vec<f64>,
    pub xmean: Vec<f64>,
}

/// Run the C7L2 simulation
pub fn run() -> Results {
    let vc: f64 = 4000.0;
    let xnt: f64 = 96.6;
    let yic: f64 = 0.0;
    let vm: f64 = 3000.0;
    let hedeg: f64 = 0.0;
    let beta: f64 = 0.8;
    let xnp: f64 = 3.0;
    let signoise: f64 = 0.001;
    let ts: f64 = 0.1;
    let num_runs: usize = 50;
    let h: f64 = 0.01;

    let normal = Normal::new(0.0, 1.0).unwrap();
    let mut rng = rand::thread_rng();

    let gfilter = 1.0 - beta * beta;
    let hfilter = (1.0 - beta).powi(2);

    let mut array_tf = Vec::new();
    let mut array_sigma = Vec::new();
    let mut array_xmean = Vec::new();

    let mut tf: f64 = 0.5;
    while tf <= 10.0 {
        let mut z = vec![0.0; num_runs];
        let mut z1 = 0.0;

        for z_item in z.iter_mut() {
            let mut y = yic;
            let mut yd = -vm * hedeg / 57.3;
            let mut t: f64 = 0.0;
            let mut s: f64 = 0.0;
            let mut xlamh = 0.0;
            let mut xlamdh = 0.0;
            let mut xnc = 0.0;

            while t <= (tf - 1e-5) {
                let y_old = y;
                let yd_old = yd;

                // First derivative evaluation
                let tgo = tf - t + 0.00001;
                let rtm = vc * tgo;
                let _xlam = y / (vc * tgo);
                let _xlamd = (rtm * yd + y * vc) / (rtm * rtm);
                let ydd = xnt - xnc;

                // Euler step
                y += h * yd;
                yd += h * ydd;
                t += h;

                // Second derivative
                let tgo = tf - t + 0.00001;
                let rtm = vc * tgo;
                let xlam = y / (vc * tgo);
                let _xlamd = (rtm * yd + y * vc) / (rtm * rtm);
                let ydd = xnt - xnc;

                // RK2 averaging
                y = 0.5 * (y_old + y + h * yd);
                yd = 0.5 * (yd_old + yd + h * ydd);

                s += h;
                if s >= (ts - 1e-5) {
                    s = 0.0;
                    // Bug fix: original used gaussc7 which doesn't exist
                    let xlamnoise = signoise * normal.sample(&mut rng);
                    let res = xlam - (xlamh + ts * xlamdh) + xlamnoise;
                    xlamh = gfilter * res + xlamh + ts * xlamdh;
                    xlamdh += hfilter * res / ts;
                    xnc = xnp * vc * xlamdh;
                }
            }

            *z_item = y;
            z1 += *z_item;
        }

        let xmean = z1 / num_runs as f64;

        // Calculate standard deviation
        let mut z1 = 0.0;
        for z_val in z.iter() {
            z1 += (*z_val - xmean).powi(2);
        }
        let sigma = if num_runs > 1 {
            (z1 / (num_runs - 1) as f64).sqrt()
        } else {
            0.0
        };

        array_tf.push(tf);
        array_sigma.push(sigma);
        array_xmean.push(xmean);

        tf += 0.5;
    }

    Results {
        tf: array_tf,
        sigma: array_sigma,
        xmean: array_xmean,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c7l2_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.tf.clone(),
        results.sigma.clone(),
        results.xmean.clone(),
    ])?;

    let plot_file = format!("{}/c7l2_sigma.png", output_dir);
    let config = PlotConfig::new("Standard deviation of miss for various flight times")
        .with_labels("Flight Time (S)", "Noise Miss Standard Deviation (Ft)")
        .with_y_range(0.0, 4.0);
    let series = vec![
        Series::new(results.tf.clone(), results.sigma.clone())
            .with_color(plotters::prelude::BLUE),
    ];
    line_plot(&plot_file, &config, &series).ok();

    let plot_file2 = format!("{}/c7l2_mean.png", output_dir);
    let config2 = PlotConfig::new("Mean of miss for various flight times")
        .with_labels("Flight Time (S)", "Mean Miss (Ft)")
        .with_y_range(0.0, 60.0);
    let series2 = vec![
        Series::new(results.tf.clone(), results.xmean.clone())
            .with_color(plotters::prelude::BLUE),
    ];
    line_plot(&plot_file2, &config2, &series2).ok();

    println!("C7L2: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c7l2_runs() {
        let results = run();
        assert!(!results.tf.is_empty());
    }
}
