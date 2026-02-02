//! Chapter 30, Lesson 1: Kalman Filter with Acceleration Estimation
//!
//! Implements an extended Kalman filter for estimating target acceleration
//! using noisy LOS angle measurements.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;
use rand::SeedableRng;
use rand_distr::{Distribution, Normal};

pub struct Results {
    pub time: Vec<f64>,
    pub xntg: Vec<f64>,
    pub xnthg: Vec<f64>,
    pub y: Vec<f64>,
    pub ystar: Vec<f64>,
}

/// Run the C30L1 simulation
pub fn run() -> Results {
    run_with_seed(None)
}

/// Run with optional seed for reproducibility
pub fn run_with_seed(seed: Option<u64>) -> Results {
    let tau: f64 = 0.5;
    let apn: i32 = 0;
    let vc: f64 = 9000.0;
    let xntreal: f64 = 96.6;
    let xntmax: f64 = 96.6;
    let w: f64 = 2.0;
    let yic: f64 = 0.0;
    let vm: f64 = 3000.0;
    let hedeg: f64 = 0.0;
    let hedegfil: f64 = 20.0;
    let xnp: f64 = 3.0;
    let sigrin: f64 = 0.001;
    let ts: f64 = 0.01;
    let tf: f64 = 10.0;

    let mut rng: rand::rngs::StdRng = match seed {
        Some(s) => rand::rngs::StdRng::seed_from_u64(s),
        None => rand::rngs::StdRng::from_entropy(),
    };
    let normal = Normal::new(0.0, 1.0).unwrap();

    let ts2 = ts * ts;
    let ts3 = ts2 * ts;
    let ts4 = ts3 * ts;
    let ts5 = ts4 * ts;
    let phin = xntmax * xntmax / tf;
    let rtm_init = vc * tf;
    let signoise = sigrin;
    let sigpos = rtm_init * signoise;
    let sign2_init = sigpos * sigpos;

    let mut p11 = sign2_init;
    let mut p12: f64 = 0.0;
    let mut p13: f64 = 0.0;
    let mut p22 = (vm * hedegfil / 57.3).powi(2);
    let mut p23: f64 = 0.0;
    let mut p33 = xntmax * xntmax;

    let mut y = yic;
    let mut yd = -vm * hedeg / 57.3;
    let mut t: f64 = 0.0;
    let h: f64 = 0.001;
    let mut s: f64 = 0.0;
    let mut yh: f64 = 0.0;
    let mut ydh: f64 = 0.0;
    let mut xnth: f64 = 0.0;
    let mut xnc: f64 = 0.0;
    let mut xnl: f64 = 0.0;

    let mut array_t = Vec::new();
    let mut array_xntg = Vec::new();
    let mut array_xnthg = Vec::new();
    let mut array_y = Vec::new();
    let mut array_ystar = Vec::new();

    while t <= tf {
        let yold = y;
        let ydold = yd;
        let xnlold = xnl;

        // First derivative
        let xnt = xntreal * (w * t).sin();
        let tgo = tf - t + 0.00001;
        let rtm = vc * tgo;
        let _xlam = y / (vc * tgo);
        let _xlamd = (rtm * yd + y * vc) / (rtm * rtm);
        let xnld = (xnc - xnl) / tau;
        let ydd = xnt - xnl;

        // Euler step
        y += h * yd;
        yd += h * ydd;
        xnl += h * xnld;
        t += h;

        // Second derivative for RK2
        let xnt = xntreal * (w * t).sin();
        let tgo = tf - t + 0.00001;
        let rtm = vc * tgo;
        let xlam = y / (vc * tgo);
        let _xlamd = (rtm * yd + y * vc) / (rtm * rtm);
        let xnld = (xnc - xnl) / tau;
        let ydd = xnt - xnl;

        // RK2 averaging
        y = 0.5 * (yold + y + h * yd);
        yd = 0.5 * (ydold + yd + h * ydd);
        xnl = 0.5 * (xnlold + xnl + h * xnld);

        s += h;
        if s >= ts - 0.0001 {
            s = 0.0;
            let tgo = tf - t + 0.000001;
            let rtm = vc * tgo;
            let signoise = sigrin;
            let sigpos = rtm * signoise;
            let sign2 = sigpos * sigpos;

            // Propagate covariance
            let m11 = p11 + ts * p12 + 0.5 * ts2 * p13
                + ts * (p12 + ts * p22 + 0.5 * ts2 * p23)
                + 0.5 * ts2 * (p13 + ts * p23 + 0.5 * ts2 * p33)
                + ts5 * phin / 20.0;
            let m12 = p12 + ts * p22 + 0.5 * ts2 * p23
                + ts * (p13 + ts * p23 + 0.5 * ts2 * p33)
                + ts4 * phin / 8.0;
            let m13 = p13 + ts * p23 + 0.5 * ts2 * p33 + phin * ts3 / 6.0;
            let m22 = p22 + ts * p23 + ts * (p23 + ts * p33) + phin * ts3 / 3.0;
            let m23 = p23 + ts * p33 + 0.5 * ts2 * phin;
            let m33 = p33 + phin * ts;

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

            // Measurement
            let xlamnoise = signoise * normal.sample(&mut rng);
            let ystar = rtm * (xlam + xlamnoise);
            let res = ystar - yh - ts * ydh - 0.5 * ts * ts * (xnth - xnl);

            // Update estimates
            yh = k1 * res + yh + ts * ydh + 0.5 * ts * ts * (xnth - xnl);
            ydh = k2 * res + ydh + ts * (xnth - xnl);
            xnth += k3 * res;

            // Guidance law
            let _xlamdh = (yh + ydh * tgo) / (vc * tgo * tgo);
            if apn == 0 {
                xnc = xnp * (yh + ydh * tgo) / (tgo * tgo);
            } else if apn == 1 {
                xnc = xnp * (yh + ydh * tgo + 0.5 * xnth * tgo * tgo) / (tgo * tgo);
            } else {
                let xs = tgo / tau;
                let top = 6.0 * xs * xs * ((-xs).exp() - 1.0 + xs);
                let bot1 = 2.0 * xs.powi(3) + 3.0 + 6.0 * xs - 6.0 * xs * xs;
                let bot2 = -12.0 * xs * (-xs).exp() - 3.0 * (-2.0 * xs).exp();
                let xnpp = top / (0.0001 + bot1 + bot2);
                let c1 = xnpp / (tgo * tgo);
                let c2 = xnpp / tgo;
                let c3 = 0.5 * xnpp;
                let c4 = -xnpp * ((-xs).exp() + xs - 1.0) / (xs * xs);
                xnc = c1 * yh + c2 * ydh + c3 * xnth + c4 * xnl;
            }

            let xntg = xnt / 32.2;
            let xnthg = xnth / 32.2;
            array_t.push(t);
            array_xntg.push(xntg);
            array_xnthg.push(xnthg);
            array_y.push(y);
            array_ystar.push(ystar);
        }
    }

    Results {
        time: array_t,
        xntg: array_xntg,
        xnthg: array_xnthg,
        y: array_y,
        ystar: array_ystar,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c30l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.xntg.clone(),
        results.xnthg.clone(),
        results.y.clone(),
        results.ystar.clone(),
    ])?;

    let plot_file = format!("{}/c30l1_accel.png", output_dir);
    let config = PlotConfig::new("Acceleration Estimate")
        .with_labels("Time (Sec)", "Acceleration (G)");

    let series = vec![
        Series::new(results.time.clone(), results.xntg.clone())
            .with_color(plotters::prelude::BLUE)
            .with_label("Actual"),
        Series::new(results.time.clone(), results.xnthg.clone())
            .with_color(plotters::prelude::RED)
            .with_label("Estimated"),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C30L1: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c30l1_runs() {
        let results = run_with_seed(Some(12345));
        assert!(!results.time.is_empty());
    }
}
