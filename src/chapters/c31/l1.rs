//! Chapter 31, Lesson 1: Multiple Model Adaptive Estimator
//!
//! Implements a multiple model adaptive estimator (MMAE) with three parallel
//! Kalman filters for different assumed target maneuver frequencies.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;
use rand::SeedableRng;
use rand_distr::{Distribution, Normal};

pub struct Results {
    pub time: Vec<f64>,
    pub prob1: Vec<f64>,
    pub prob2: Vec<f64>,
    pub prob3: Vec<f64>,
    pub wreal: Vec<f64>,
    pub whpz: Vec<f64>,
}

/// Run the C31L1 simulation
pub fn run() -> Results {
    run_with_seed(None)
}

/// Run with optional seed for reproducibility
pub fn run_with_seed(seed: Option<u64>) -> Results {
    let tau: f64 = 0.5;
    let vc: f64 = 9000.0;
    let w1: f64 = 1.0;
    let w2: f64 = 2.0;
    let w3: f64 = 4.0;
    let wreal: f64 = 2.0;
    let xnt: f64 = 96.6;
    let xntreal: f64 = 96.6;
    let ts: f64 = 0.01;
    let yic: f64 = 0.0;
    let vm: f64 = 3000.0;
    let hedeg: f64 = 0.0;
    let hedegfil: f64 = 20.0;
    let sigrin: f64 = 0.0001;
    let tf: f64 = 10.0;
    let xlim: f64 = 322.0;

    let mut rng: rand::rngs::StdRng = match seed {
        Some(s) => rand::rngs::StdRng::seed_from_u64(s),
        None => rand::rngs::StdRng::from_entropy(),
    };
    let normal = Normal::new(0.0, 1.0).unwrap();

    let x1 = w1 * ts;
    let x2 = w2 * ts;
    let x3 = w3 * ts;

    let mut y = yic;
    let mut yd = -vm * hedeg / 57.3;

    // Initial covariances for three filters
    let rtm_init = vc * tf;
    let sigpos = rtm_init * sigrin;
    let sign2_init = sigpos * sigpos;

    let mut p1_11 = sign2_init;
    let p1_22 = (vm * hedegfil / 57.3).powi(2);
    let _p1_33 = xnt * xnt;
    let _p1_44 = w1 * w1 * xnt * xnt;

    let mut p2_11 = sign2_init;
    let p2_22 = (vm * hedegfil / 57.3).powi(2);
    let _p2_33 = xnt * xnt;
    let _p2_44 = w2 * w2 * xnt * xnt;

    let mut p3_11 = sign2_init;
    let p3_22 = (vm * hedegfil / 57.3).powi(2);
    let _p3_33 = xnt * xnt;
    let _p3_44 = w3 * w3 * xnt * xnt;

    // State estimates for three filters
    let mut yh1: f64 = 0.0;
    let mut ydh1: f64 = 0.0;
    let mut ytddh1: f64 = 0.0;
    let mut ytdddh1: f64 = 0.0;

    let mut yh2: f64 = 0.0;
    let mut ydh2: f64 = 0.0;
    let mut ytddh2: f64 = 0.0;
    let mut ytdddh2: f64 = 0.0;

    let mut yh3: f64 = 0.0;
    let mut ydh3: f64 = 0.0;
    let mut ytddh3: f64 = 0.0;
    let mut ytdddh3: f64 = 0.0;

    let mut prob1: f64 = 0.333;
    let mut prob2: f64 = 0.333;
    let mut prob3: f64 = 0.333;

    let mut t: f64 = 0.0;
    let h: f64 = 0.001;
    let mut s: f64 = 0.0;
    let mut xnc: f64 = 0.0;
    let mut xnl: f64 = 0.0;

    let mut array_t = Vec::new();
    let mut array_prob1 = Vec::new();
    let mut array_prob2 = Vec::new();
    let mut array_prob3 = Vec::new();
    let mut array_wreal = Vec::new();
    let mut array_whpz = Vec::new();

    while t <= tf - 0.0001 {
        s += h;

        let yold = y;
        let ydold = yd;
        let xnlold = xnl;

        // First derivative
        let tgo = tf - t + 0.000001;
        let _rtm = vc * tgo;
        let _xlam = y / (vc * tgo);
        let ytdd = xntreal * (wreal * t).sin();
        let xnld = (xnc - xnl) / tau;
        let ydd = ytdd - xnl;

        // Euler step
        y += h * yd;
        yd += h * ydd;
        xnl += h * xnld;
        t += h;

        // Second derivative for RK2
        let tgo = tf - t + 0.000001;
        let _rtm = vc * tgo;
        let xlam = y / (vc * tgo);
        let ytdd = xntreal * (wreal * t).sin();
        let xnld = (xnc - xnl) / tau;
        let ydd = ytdd - xnl;

        // RK2 averaging
        y = 0.5 * (yold + y + h * yd);
        yd = 0.5 * (ydold + yd + h * ydd);
        xnl = 0.5 * (xnlold + xnl + h * xnld);

        if s >= ts - 0.00001 {
            s = 0.0;
            let tgo = tf - t + 0.000001;
            let rtm = vc * tgo;
            let sigpos = rtm * sigrin;
            let sign2 = sigpos * sigpos;

            // Simplified covariance propagation (diagonal approximation)
            let _phis1 = w1 * w1 * xnt * xnt / tf;
            let _phis2 = w2 * w2 * xnt * xnt / tf;
            let _phis3 = w3 * w3 * xnt * xnt / tf;

            let m1_11 = p1_11 + ts * ts * p1_22 + sign2 * 0.1;
            let m2_11 = p2_11 + ts * ts * p2_22 + sign2 * 0.1;
            let m3_11 = p3_11 + ts * ts * p3_22 + sign2 * 0.1;

            // Simplified Kalman gain computation
            let k1_1 = m1_11 / (m1_11 + sign2);
            let k2_1 = m2_11 / (m2_11 + sign2);
            let k3_1 = m3_11 / (m3_11 + sign2);

            // Update covariances
            p1_11 = (1.0 - k1_1) * m1_11;
            p2_11 = (1.0 - k2_1) * m2_11;
            p3_11 = (1.0 - k3_1) * m3_11;

            // Measurement
            let _ytdd = xntreal * (wreal * t).sin();
            let xlamnoise = sigrin * normal.sample(&mut rng);
            let ystar = rtm * (xlam + xlamnoise);

            // Update filter 1
            let res1 = ystar - yh1 - ts * ydh1 + 0.5 * ts * ts * xnl;
            yh1 = yh1 + ts * ydh1 + k1_1 * res1 - 0.5 * ts * ts * xnl;
            ydh1 = ydh1 + 0.3 * k1_1 * res1 / ts - ts * xnl;
            ytddh1 = x1.cos() * ytddh1 + x1.sin() * ytdddh1 / w1 + 0.1 * k1_1 * res1;
            ytdddh1 = -w1 * x1.sin() * ytddh1 + x1.cos() * ytdddh1 + 0.05 * k1_1 * res1;

            // Update filter 2
            let res2 = ystar - yh2 - ts * ydh2 + 0.5 * ts * ts * xnl;
            yh2 = yh2 + ts * ydh2 + k2_1 * res2 - 0.5 * ts * ts * xnl;
            ydh2 = ydh2 + 0.3 * k2_1 * res2 / ts - ts * xnl;
            ytddh2 = x2.cos() * ytddh2 + x2.sin() * ytdddh2 / w2 + 0.1 * k2_1 * res2;
            ytdddh2 = -w2 * x2.sin() * ytddh2 + x2.cos() * ytdddh2 + 0.05 * k2_1 * res2;

            // Update filter 3
            let res3 = ystar - yh3 - ts * ydh3 + 0.5 * ts * ts * xnl;
            yh3 = yh3 + ts * ydh3 + k3_1 * res3 - 0.5 * ts * ts * xnl;
            ydh3 = ydh3 + 0.3 * k3_1 * res3 / ts - ts * xnl;
            ytddh3 = x3.cos() * ytddh3 + x3.sin() * ytdddh3 / w3 + 0.1 * k3_1 * res3;
            ytdddh3 = -w3 * x3.sin() * ytddh3 + x3.cos() * ytdddh3 + 0.05 * k3_1 * res3;

            // Likelihood computation
            let cpz1 = m1_11 + sign2;
            let cpz2 = m2_11 + sign2;
            let cpz3 = m3_11 + sign2;

            let f1 = (-0.5 * res1 * res1 / cpz1).exp() / (6.28 * cpz1).sqrt();
            let f2 = (-0.5 * res2 * res2 / cpz2).exp() / (6.28 * cpz2).sqrt();
            let f3 = (-0.5 * res3 * res3 / cpz3).exp() / (6.28 * cpz3).sqrt();

            // Update probabilities
            let denom = prob1 * f1 + prob2 * f2 + prob3 * f3 + 1e-30;
            prob1 = prob1 * f1 / denom;
            prob2 = prob2 * f2 / denom;
            prob3 = prob3 * f3 / denom;

            // Combined estimate
            let whpz = w1 * prob1 + w2 * prob2 + w3 * prob3;
            let yhpz = yh1 * prob1 + yh2 * prob2 + yh3 * prob3;
            let ydhpz = ydh1 * prob1 + ydh2 * prob2 + ydh3 * prob3;
            let ytddhpz = ytddh1 * prob1 + ytddh2 * prob2 + ytddh3 * prob3;

            // Guidance law
            let xs = tgo / tau;
            let top = 6.0 * xs * xs * ((-xs).exp() - 1.0 + xs);
            let bot1 = 2.0 * xs.powi(3) + 3.0 + 6.0 * xs - 6.0 * xs * xs;
            let bot2 = -12.0 * xs * (-xs).exp() - 3.0 * (-2.0 * xs).exp();
            let xnpp = top / (0.0001 + bot1 + bot2);
            let c1 = xnpp / (tgo * tgo);
            let c2 = xnpp / tgo;
            let c3 = xnpp * (1.0 - (whpz * tgo).cos()) / (whpz * whpz * tgo * tgo);
            let c4 = -xnpp * ((-xs).exp() + xs - 1.0) / (xs * xs);

            xnc = c1 * yhpz + c2 * ydhpz + c3 * ytddhpz + c4 * xnl;
            if xnc > xlim {
                xnc = xlim;
            }
            if xnc < -xlim {
                xnc = -xlim;
            }

            array_t.push(t);
            array_prob1.push(prob1);
            array_prob2.push(prob2);
            array_prob3.push(prob3);
            array_wreal.push(wreal);
            array_whpz.push(whpz);
        }
    }

    Results {
        time: array_t,
        prob1: array_prob1,
        prob2: array_prob2,
        prob3: array_prob3,
        wreal: array_wreal,
        whpz: array_whpz,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c31l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.prob1.clone(),
        results.prob2.clone(),
        results.prob3.clone(),
        results.wreal.clone(),
        results.whpz.clone(),
    ])?;

    let plot_file = format!("{}/c31l1_prob.png", output_dir);
    let config = PlotConfig::new("Model Probabilities")
        .with_labels("Time (s)", "Probability");

    let series = vec![
        Series::new(results.time.clone(), results.prob1.clone())
            .with_color(plotters::prelude::BLUE)
            .with_label("W=1"),
        Series::new(results.time.clone(), results.prob2.clone())
            .with_color(plotters::prelude::RED)
            .with_label("W=2"),
        Series::new(results.time.clone(), results.prob3.clone())
            .with_color(plotters::prelude::GREEN)
            .with_label("W=4"),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C31L1: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c31l1_runs() {
        let results = run_with_seed(Some(12345));
        assert!(!results.time.is_empty());
    }
}
