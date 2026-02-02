//! Chapter 6, Lesson 3: Fifth-Order Binomial Guidance with Rate Feedback
//!
//! Simulates missile-target engagement with fifth-order binomial guidance
//! and rate feedback compensation.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub tf: Vec<f64>,
    pub y: Vec<f64>,
}

/// Run the C6L3 simulation
pub fn run() -> Results {
    let vc: f64 = 4000.0;
    let xnt: f64 = 32.2;
    let yic: f64 = 0.0;
    let vm: f64 = 3000.0;
    let hedeg: f64 = 0.0;
    let tau: f64 = 0.5;
    let xnp: f64 = 3.0;
    let ta: f64 = 0.0;
    let r: f64 = -0.01;
    let h: f64 = 0.01;

    let mut array_tf = Vec::new();
    let mut array_y = Vec::new();

    let mut tf: f64 = 0.1;
    while tf <= 10.0 {
        let mut y = yic;
        let mut yd = -vm * hedeg / 57.3;
        let mut xnl = 0.0;
        let mut elamdh = 0.0;
        let mut x4 = 0.0;
        let mut x5 = 0.0;
        let mut th = 0.0;
        let mut thh = 0.0;
        let mut t: f64 = 0.0;

        while t <= (tf - 1e-5) {
            let y_old = y;
            let yd_old = yd;
            let xnl_old = xnl;
            let elamdh_old = elamdh;
            let x4_old = x4;
            let x5_old = x5;
            let th_old = th;
            let thh_old = thh;

            // First derivative evaluation
            let tgo = tf - t + 0.00001;
            let xlam = y / (vc * tgo);
            let eps = xlam - th - thh + r * thh;
            let dd = 5.0 * eps / tau;
            let elamdhd = 5.0 * (dd - elamdh) / tau;
            let xnc = xnp * vc * elamdh;
            let x4d = 5.0 * (xnc - x4) / tau;
            let x5d = 5.0 * (x4 - x5) / tau;
            let xnld = 5.0 * (x5 - xnl) / tau;
            let thd = xnl / vm + ta * xnld / vm;
            let thhd = dd - thd;
            let ydd = xnt - xnl;

            // Euler step
            y += h * yd;
            yd += h * ydd;
            xnl += h * xnld;
            elamdh += h * elamdhd;
            x4 += h * x4d;
            x5 += h * x5d;
            th += h * thd;
            thh += h * thhd;
            t += h;

            // Second derivative evaluation
            let tgo = tf - t + 0.00001;
            let xlam = y / (vc * tgo);
            let eps = xlam - th - thh + r * thh;
            let dd = 5.0 * eps / tau;
            let elamdhd = 5.0 * (dd - elamdh) / tau;
            let xnc = xnp * vc * elamdh;
            let x4d = 5.0 * (xnc - x4) / tau;
            let x5d = 5.0 * (x4 - x5) / tau;
            let xnld = 5.0 * (x5 - xnl) / tau;
            let thd = xnl / vm + ta * xnld / vm;
            let thhd = dd - thd;
            let ydd = xnt - xnl;

            // RK2 averaging
            y = 0.5 * (y_old + y + h * yd);
            yd = 0.5 * (yd_old + yd + h * ydd);
            xnl = 0.5 * (xnl_old + xnl + h * xnld);
            elamdh = 0.5 * (elamdh_old + elamdh + h * elamdhd);
            x4 = 0.5 * (x4_old + x4 + h * x4d);
            x5 = 0.5 * (x5_old + x5 + h * x5d);
            th = 0.5 * (th_old + th + h * thd);
            thh = 0.5 * (thh_old + thh + h * thhd);
        }

        array_tf.push(tf);
        array_y.push(y);

        tf += 0.1;
    }

    Results {
        tf: array_tf,
        y: array_y,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c6l3_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.tf.clone(),
        results.y.clone(),
    ])?;

    let plot_file = format!("{}/c6l3_plot.png", output_dir);
    let config = PlotConfig::new("Fifth-Order Binomial Guidance Miss")
        .with_labels("Flight Time (Sec)", "Miss (Ft)");
    let series = vec![
        Series::new(results.tf.clone(), results.y.clone())
            .with_color(plotters::prelude::BLUE),
    ];
    line_plot(&plot_file, &config, &series).ok();

    println!("C6L3: Simulation finished");
    println!("  Data saved to: {}", data_file);
    println!("  Plot saved to: {}", plot_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c6l3_runs() {
        let results = run();
        assert!(!results.tf.is_empty());
    }
}
