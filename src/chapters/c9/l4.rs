//! Chapter 9, Lesson 4: Kalman Filter with Variable Target Maneuver
//!
//! Kalman filter simulation with step, sinusoidal, or square wave target maneuver.
//! Note: This is a deterministic version (XLAMNOISE=0) for verification purposes.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub xntg: Vec<f64>,
    pub xnthg: Vec<f64>,
}

/// Run the C9L4 simulation
pub fn run() -> Results {
    let iconstant: i32 = 1;  // 1=step, 2=sin, 3=square
    let w: f64 = 2.0;
    let tstart: f64 = 3.0;
    let vc: f64 = 4000.0;
    let xntic: f64 = 96.6;
    let yic: f64 = 0.0;
    let vm: f64 = 3000.0;
    let hedeg: f64 = 0.0;
    let hedegfil: f64 = 20.0;
    let xnp: f64 = 3.0;
    let sigrin: f64 = 0.001;
    let ts: f64 = 0.1;
    let apn: f64 = 0.0;
    let tf: f64 = 10.0;
    let h: f64 = 0.001;

    let ts2 = ts * ts;
    let ts3 = ts2 * ts;
    let ts4 = ts3 * ts;
    let ts5 = ts4 * ts;
    let phin = xntic * xntic / tf;
    let rtm_init = vc * tf;
    let signoise = sigrin;
    let sigpos = rtm_init * signoise;
    let sign2_init = sigpos * sigpos;

    let mut p11 = sign2_init;
    let mut p12 = 0.0;
    let mut p13 = 0.0;
    let mut p22 = (vm * hedegfil / 57.3).powi(2);
    let mut p23 = 0.0;
    let mut p33 = xntic * xntic;

    let mut y = yic;
    let mut yd = -vm * hedeg / 57.3;
    let mut t: f64 = 0.0;
    let mut s: f64 = 0.0;

    let mut yh = 0.0;
    let mut ydh = 0.0;
    let mut xnth = 0.0;
    let mut xnc = 0.0;

    let mut array_t = Vec::new();
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

        // Target acceleration based on maneuver type
        let xntc = match iconstant {
            1 => if t < tstart { 0.0 } else { xntic },
            2 => if t < tstart { 0.0 } else { xntic * (w * t).sin() },
            _ => if t < tstart { 0.0 } else { xntic * (w * (t - tstart)).sin().signum() },
        };

        let ydd = xntc - xnc;

        // Euler step
        y += h * yd;
        yd += h * ydd;
        t += h;

        // Second derivative
        let tgo = tf - t + 0.00001;
        let rtm = vc * tgo;
        let xlam = y / (vc * tgo);
        let _xlamd = (rtm * yd + y * vc) / (rtm * rtm);

        let xntc = match iconstant {
            1 => if t < tstart { 0.0 } else { xntic },
            2 => if t < tstart { 0.0 } else { xntic * (w * t).sin() },
            _ => if t < tstart { 0.0 } else { xntic * (w * (t - tstart)).sin().signum() },
        };

        let ydd = xntc - xnc;

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

            let k1 = m11 / (m11 + sign2);
            let k2 = m12 / (m11 + sign2);
            let k3 = m13 / (m11 + sign2);

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

            yh = k1 * res + yh + ts * ydh + 0.5 * ts * ts * (xnth - xnc);
            ydh = k2 * res + ydh + ts * (xnth - xnc);
            xnth += k3 * res;

            let xlamdh = (yh + ydh * tgo) / (vc * tgo * tgo);
            xnc = xnp * vc * xlamdh + apn * 0.5 * xnp * xnth;

            array_t.push(t);
            array_xntg.push(xntc / 32.2);
            array_xnthg.push(xnth / 32.2);
        }
    }

    Results {
        time: array_t,
        xntg: array_xntg,
        xnthg: array_xnthg,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c9l4_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.xntg.clone(),
        results.xnthg.clone(),
    ])?;

    let plot_file = format!("{}/c9l4_plot.png", output_dir);
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

    println!("C9L4: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c9l4_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
