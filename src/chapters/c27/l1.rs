//! Chapter 27, Lesson 1: Fading Memory Filter Miss Analysis
//!
//! Computes standard miss distance for various flight times using
//! a fading memory filter for LOS rate estimation.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub tf: Vec<f64>,
    pub y_miss: Vec<f64>,
}

/// Run the C27L1 simulation
pub fn run() -> Results {
    let vc: f64 = 4000.0;
    let tap: f64 = 0.2;
    let xnt: f64 = 96.6;
    let beta: f64 = 0.8;
    let xnp: f64 = 3.0;
    let ts: f64 = 0.1;
    let ts2: f64 = 0.02;

    let mut array_tf = Vec::new();
    let mut array_y = Vec::new();

    let mut tf_val = 0.1;
    while tf_val <= 10.0 {
        let mut y: f64 = 0.0;
        let mut yd: f64 = 0.0;
        let mut t: f64 = 0.0;
        let h: f64 = 0.001;
        let mut s: f64 = 0.0;
        let mut s2: f64 = 0.0;

        let gfilter = 1.0 - beta * beta;
        let hfilter = (1.0 - beta).powi(2);

        let mut xlamhold: f64 = 0.0;
        let mut xlamdhold: f64 = 0.0;
        let mut y1old: f64 = 0.0;
        let mut y2old: f64 = 0.0;
        let mut y3old: f64 = 0.0;
        let mut y4old: f64 = 0.0;
        let mut y5old: f64 = 0.0;
        let mut xnc: f64 = 0.0;
        let mut xnl: f64 = 0.0;
        let mut pz: f64 = 0.0;
        let mut xlam: f64;

        while t <= tf_val - 1e-5 {
            let yold = y;
            let ydold = yd;
            let xnlold = xnl;

            // First derivative evaluation
            let tgo = tf_val - t + 0.00001;
            let rtm = vc * tgo;
            let _xlam = y / (vc * tgo);
            let _xlamd = (rtm * yd + y * vc) / (rtm * rtm);
            let xnld = (xnc - xnl) / tap;
            let ydd = xnt - xnl;

            // Euler step
            y += h * yd;
            yd += h * ydd;
            xnl += h * xnld;
            t += h;

            // Second derivative for RK2
            let tgo = tf_val - t + 0.00001;
            let rtm = vc * tgo;
            xlam = y / (vc * tgo);
            let _xlamd = (rtm * yd + y * vc) / (rtm * rtm);
            let xnld = (xnc - xnl) / tap;
            let ydd = xnt - xnl;

            // RK2 averaging
            y = 0.5 * (yold + y + h * yd);
            yd = 0.5 * (ydold + yd + h * ydd);
            xnl = 0.5 * (xnlold + xnl + h * xnld);

            s += h;
            s2 += h;

            if s2 >= ts2 - 1e-5 {
                s2 = 0.0;
                let y1new = xlam;
                let y2new = y1old;
                let y3new = y2old;
                let y4new = y3old;
                let y5new = y4old;
                pz = 0.2 * (y5old + y5new + y4new + y3new + y2new + xlam);
                y1old = y1new;
                y2old = y2new;
                y3old = y3new;
                y4old = y4new;
                y5old = y5new;
            }

            if s >= ts - 1e-5 {
                s = 0.0;
                let res = pz - (xlamhold + ts * xlamdhold);
                let xlamhnew = gfilter * res + xlamhold + ts * xlamdhold;
                let xlamdhnew = hfilter * res / ts + xlamdhold;
                xnc = xnp * vc * xlamdhnew;
                xlamhold = xlamhnew;
                xlamdhold = xlamdhnew;
            }
        }

        array_tf.push(tf_val);
        array_y.push(y);

        tf_val += 0.1;
    }

    Results {
        tf: array_tf,
        y_miss: array_y,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c27l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.tf.clone(),
        results.y_miss.clone(),
    ])?;

    let plot_file = format!("{}/c27l1_miss.png", output_dir);
    let config = PlotConfig::new("Standard Miss for Various Flight Times")
        .with_labels("Flight Time (S)", "Miss (Ft)");

    let series = vec![
        Series::new(results.tf.clone(), results.y_miss.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C27L1: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c27l1_runs() {
        let results = run();
        assert!(!results.tf.is_empty());
    }
}
