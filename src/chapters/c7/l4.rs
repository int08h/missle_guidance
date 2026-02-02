//! Chapter 7, Lesson 4: Third-Order Digital Filter
//!
//! Simulation with third-order fading memory filter.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;
use rand_distr::{Distribution, Normal};

pub struct Results {
    pub time: Vec<f64>,
    pub y: Vec<f64>,
    pub xncg: Vec<f64>,
    pub xlamd: Vec<f64>,
    pub xlamdh: Vec<f64>,
    pub xntg: Vec<f64>,
    pub xnthg: Vec<f64>,
}

/// Run the C7L4 simulation
pub fn run() -> Results {
    let vc: f64 = 4000.0;
    let xnt: f64 = 96.6;
    let yic: f64 = 0.0;
    let vm: f64 = 3000.0;
    let hedeg: f64 = 0.0;
    let beta: f64 = 0.8;
    let xnp: f64 = 3.0;
    let signoise: f64 = 0.001;
    let tf: f64 = 10.0;
    let ts: f64 = 0.1;
    let h: f64 = 0.01;

    let normal = Normal::new(0.0, 1.0).unwrap();
    let mut rng = rand::thread_rng();

    // Third-order filter coefficients
    let gfilter = 1.0 - beta.powi(3);
    let hfilter = 1.5 * (1.0 - beta).powi(2) * (1.0 + beta);
    let kfilter = 0.5 * (1.0 - beta).powi(3);

    let mut y = yic;
    let mut yd = -vm * hedeg / 57.3;
    let mut t: f64 = 0.0;
    let mut s: f64 = 0.0;

    let mut yh = 0.0;
    let mut ydh = 0.0;
    let mut xnth = 0.0;
    let mut xnc = 0.0;

    let mut array_t = Vec::new();
    let mut array_y = Vec::new();
    let mut array_xncg = Vec::new();
    let mut array_xlamd = Vec::new();
    let mut array_xlamdh = Vec::new();
    let mut array_xntg = Vec::new();
    let mut array_xnthg = Vec::new();

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
        let xlamd_new = (rtm * yd + y * vc) / (rtm * rtm);
        let ydd = xnt - xnc;

        // RK2 averaging
        y = 0.5 * (y_old + y + h * yd);
        yd = 0.5 * (yd_old + yd + h * ydd);

        s += h;
        if s >= (ts - 1e-5) {
            s = 0.0;
            let xlamnoise = signoise * normal.sample(&mut rng);
            let ystar = rtm * (xlam + xlamnoise);
            let res = ystar - yh - ts * ydh - 0.5 * ts * ts * (xnth - xnc);
            yh = gfilter * res + yh + ts * ydh + 0.5 * ts * ts * (xnth - xnc);
            ydh = hfilter * res / ts + ydh + ts * (xnth - xnc);
            xnth += 2.0 * kfilter * res / (ts * ts);
            let xlamdh = (yh + ydh * tgo) / (vc * tgo * tgo);
            xnc = xnp * vc * xlamdh;

            array_t.push(t);
            array_y.push(y);
            array_xncg.push(xnc / 32.2);
            array_xlamd.push(xlamd_new);
            array_xlamdh.push(xlamdh);
            array_xntg.push(xnt / 32.2);
            array_xnthg.push(xnth / 32.2);
        }
    }

    Results {
        time: array_t,
        y: array_y,
        xncg: array_xncg,
        xlamd: array_xlamd,
        xlamdh: array_xlamdh,
        xntg: array_xntg,
        xnthg: array_xnthg,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c7l4_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.y.clone(),
        results.xncg.clone(),
        results.xlamd.clone(),
        results.xlamdh.clone(),
        results.xntg.clone(),
        results.xnthg.clone(),
    ])?;

    let plot_file = format!("{}/c7l4_los_rate.png", output_dir);
    let config = PlotConfig::new("Line of Sight Rate")
        .with_labels("Time (S)", "Line of Sight Rate (Rad/S)")
        .with_y_range(0.0, 0.05);
    let series = vec![
        Series::new(results.time.clone(), results.xlamd.clone())
            .with_label("Actual")
            .with_color(plotters::prelude::BLUE),
        Series::new(results.time.clone(), results.xlamdh.clone())
            .with_label("Estimated")
            .with_color(plotters::prelude::RED),
    ];
    line_plot(&plot_file, &config, &series).ok();

    let plot_file2 = format!("{}/c7l4_accel.png", output_dir);
    let config2 = PlotConfig::new("Acceleration")
        .with_labels("Time (S)", "Acceleration (G)");
    let series2 = vec![
        Series::new(results.time.clone(), results.xntg.clone())
            .with_label("Target")
            .with_color(plotters::prelude::BLUE),
        Series::new(results.time.clone(), results.xnthg.clone())
            .with_label("Estimated")
            .with_color(plotters::prelude::RED),
    ];
    line_plot(&plot_file2, &config2, &series2).ok();

    println!("C7L4: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c7l4_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
