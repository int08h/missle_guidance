//! Chapter 41, Lesson 1: Kalman Filter Miss Estimation
//!
//! Covariance-based estimation of RMS miss distance.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub form: Vec<f64>,
    pub sp11: Vec<f64>,
}

/// Run the C41L1 simulation
pub fn run() -> Results {
    let vc: f64 = 5.0 * 3280.0;
    let xntic: f64 = 161.0;
    let vm: f64 = 3000.0;
    let hedeg: f64 = 20.0;
    let signoise: f64 = 10.0;
    let ts: f64 = 0.01;

    let ts2 = ts * ts;
    let ts3 = ts2 * ts;
    let ts4 = ts3 * ts;
    let ts5 = ts4 * ts;

    let phin = signoise * signoise * ts;

    let mut array_t = Vec::new();
    let mut array_form = Vec::new();
    let mut array_sp11 = Vec::new();

    let mut tf = 0.2;
    while tf <= 10.0 {
        let phis = xntic * xntic / tf;
        let sigpos = signoise;
        let sign2 = sigpos * sigpos;

        let mut p11 = sign2;
        let mut p12: f64 = 0.0;
        let mut p13: f64 = 0.0;
        let mut p22 = (vm * hedeg / 57.3).powi(2);
        let mut p23: f64 = 0.0;
        let mut p33 = xntic * xntic;

        let mut t = ts;
        while t <= tf {
            let _tgo = tf - t + 0.000001;
            let _rtm = vc * _tgo;
            let sigpos = signoise;
            let sign2 = sigpos * sigpos;

            // Propagate covariance
            let m11 = p11 + ts * p12 + 0.5 * ts2 * p13
                + ts * (p12 + ts * p22 + 0.5 * ts2 * p23)
                + 0.5 * ts2 * (p13 + ts * p23 + 0.5 * ts2 * p33)
                + ts5 * phis / 20.0;
            let m12 = p12 + ts * p22 + 0.5 * ts2 * p23
                + ts * (p13 + ts * p23 + 0.5 * ts2 * p33)
                + ts4 * phis / 8.0;
            let m13 = p13 + ts * p23 + 0.5 * ts2 * p33 + phis * ts3 / 6.0;
            let m22 = p22 + ts * p23 + ts * (p23 + ts * p33) + phis * ts3 / 3.0;
            let m23 = p23 + ts * p33 + 0.5 * ts2 * phis;
            let m33 = p33 + phis * ts;

            // Kalman gains
            let k1 = m11 / (m11 + sign2);
            let k2 = m12 / (m11 + sign2);
            let k3 = m13 / (m11 + sign2);

            // Update covariance
            p11 = (1.0 - k1) * m11;
            p12 = (1.0 - k1) * m12;
            p13 = (1.0 - k1) * m13;
            p22 = -k2 * m12 + m22;
            p23 = -k2 * m13 + m23;
            p33 = -k3 * m13 + m33;

            t += ts;
        }

        let sp11 = p11.sqrt();
        let form = (2.0 * phis.powf(0.16667) * phin.powf(0.83333)).sqrt();

        array_t.push(t);
        array_form.push(form);
        array_sp11.push(sp11);

        tf += 0.2;
    }

    Results {
        time: array_t,
        form: array_form,
        sp11: array_sp11,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c41l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.form.clone(),
        results.sp11.clone(),
    ])?;

    let plot_file = format!("{}/c41l1_miss.png", output_dir);
    let config = PlotConfig::new("RMS Miss Estimation")
        .with_labels("Flight Time (S)", "RMS Miss (Ft)");

    let series = vec![
        Series::new(results.time.clone(), results.form.clone())
            .with_color(plotters::prelude::BLUE)
            .with_label("Formula"),
        Series::new(results.time.clone(), results.sp11.clone())
            .with_color(plotters::prelude::RED)
            .with_label("Kalman"),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C41L1: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c41l1_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
