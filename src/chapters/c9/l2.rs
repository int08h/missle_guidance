//! Chapter 9, Lesson 2: Kalman Filter with Polynomial Dynamics
//!
//! Kalman filter for missile guidance with polynomial target acceleration model.
//! Note: This is a deterministic version (XLAMNOISE=0) for verification purposes.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub y: Vec<f64>,
    pub xncg: Vec<f64>,
    pub xntg: Vec<f64>,
    pub xnthg: Vec<f64>,
    pub errntg: Vec<f64>,
    pub sp33g: Vec<f64>,
    pub sp33pg: Vec<f64>,
}

/// Run the C9L2 simulation
pub fn run() -> Results {
    let vc: f64 = 4000.0;
    let xnt: f64 = 96.6;
    let yic: f64 = 0.0;
    let vm: f64 = 3000.0;
    let hedeg: f64 = 0.0;
    let hedegfil: f64 = 20.0;
    let xnp: f64 = 3.0;
    let sigrin: f64 = 0.001;
    let ts: f64 = 0.1;
    let apn: f64 = 0.0;
    let tf: f64 = 10.0;
    let h: f64 = 0.01;

    let ts2 = ts * ts;
    let ts3 = ts2 * ts;
    let ts4 = ts3 * ts;
    let ts5 = ts4 * ts;
    let phin = xnt * xnt / tf;
    let rtm_init = vc * tf;
    let signoise = sigrin;
    let sigpos = rtm_init * signoise;
    let sign2_init = sigpos * sigpos;

    // Covariance matrix initialization
    let mut p11 = sign2_init;
    let mut p12 = 0.0;
    let mut p13 = 0.0;
    let mut p22 = (vm * hedegfil / 57.3).powi(2);
    let mut p23 = 0.0;
    let mut p33 = xnt * xnt;

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
    let mut array_xntg = Vec::new();
    let mut array_xnthg = Vec::new();
    let mut array_errntg = Vec::new();
    let mut array_sp33g = Vec::new();
    let mut array_sp33pg = Vec::new();

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
            let tgo = tf - t + 0.000001;
            let rtm = vc * tgo;
            let signoise = sigrin;
            let sigpos = rtm * signoise;
            let sign2 = sigpos * sigpos;

            // Covariance propagation
            let m11 = p11 + ts * p12 + 0.5 * ts2 * p13 + ts * (p12 + ts * p22 + 0.5 * ts2 * p23);
            let m11 = m11 + 0.5 * ts2 * (p13 + ts * p23 + 0.5 * ts2 * p33) + ts5 * phin / 20.0;
            let m12 = p12 + ts * p22 + 0.5 * ts2 * p23 + ts * (p13 + ts * p23 + 0.5 * ts2 * p33) + ts4 * phin / 8.0;
            let m13 = p13 + ts * p23 + 0.5 * ts2 * p33 + phin * ts3 / 6.0;
            let m22 = p22 + ts * p23 + ts * (p23 + ts * p33) + phin * ts3 / 3.0;
            let m23 = p23 + ts * p33 + 0.5 * ts2 * phin;
            let m33 = p33 + phin * ts;

            // Kalman gains
            let k1 = m11 / (m11 + sign2);
            let k2 = m12 / (m11 + sign2);
            let k3 = m13 / (m11 + sign2);

            // Covariance update
            p11 = (1.0 - k1) * m11;
            p12 = (1.0 - k1) * m12;
            p13 = (1.0 - k1) * m13;
            p22 = -k2 * m12 + m22;
            p23 = -k2 * m13 + m23;
            p33 = -k3 * m13 + m33;

            // Measurement (no noise for verification - matches C9L5 pattern)
            let xlamnoise = 0.0;
            let ystar = rtm * (xlam + xlamnoise);
            let res = ystar - yh - ts * ydh - 0.5 * ts * ts * (xnth - xnc);

            // State update
            yh = k1 * res + yh + ts * ydh + 0.5 * ts * ts * (xnth - xnc);
            ydh = k2 * res + ydh + ts * (xnth - xnc);
            xnth += k3 * res;

            let xlamdh = (yh + ydh * tgo) / (vc * tgo * tgo);
            xnc = xnp * vc * xlamdh + apn * 0.5 * xnp * xnth;

            let errnt = xnt - xnth;
            let sp33 = p33.sqrt();

            array_t.push(t);
            array_y.push(y);
            array_xncg.push(xnc / 32.2);
            array_xntg.push(xnt / 32.2);
            array_xnthg.push(xnth / 32.2);
            array_errntg.push(errnt / 32.2);
            array_sp33g.push(sp33 / 32.2);
            array_sp33pg.push(-sp33 / 32.2);
        }
    }

    Results {
        time: array_t,
        y: array_y,
        xncg: array_xncg,
        xntg: array_xntg,
        xnthg: array_xnthg,
        errntg: array_errntg,
        sp33g: array_sp33g,
        sp33pg: array_sp33pg,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c9l2_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.y.clone(),
        results.xncg.clone(),
        results.xntg.clone(),
        results.xnthg.clone(),
        results.errntg.clone(),
        results.sp33g.clone(),
        results.sp33pg.clone(),
    ])?;

    let plot_file = format!("{}/c9l2_accel.png", output_dir);
    let config = PlotConfig::new("Acceleration Estimation")
        .with_labels("Time (S)", "Acceleration (G)");
    let series = vec![
        Series::new(results.time.clone(), results.xntg.clone())
            .with_label("True")
            .with_color(plotters::prelude::BLUE),
        Series::new(results.time.clone(), results.xnthg.clone())
            .with_label("Estimated")
            .with_color(plotters::prelude::RED),
    ];
    line_plot(&plot_file, &config, &series).ok();

    println!("C9L2: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c9l2_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
