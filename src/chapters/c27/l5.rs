//! Chapter 27, Lesson 5: Miss for Various TGO with Blind Range
//!
//! Similar to L3 but includes blind range consideration (TBLIND).
//! When TGO < TBLIND, measurement noise becomes very large (blind).

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub tgos: Vec<f64>,
    pub y: Vec<f64>,
}

/// Run the C27L5 simulation
pub fn run() -> Results {
    let qzero: i32 = 0;
    let tblind: f64 = 0.5;
    let vm: f64 = 3000.0;
    let vc: f64 = 4000.0;
    let xnt: f64 = 96.6;
    let yic: f64 = 0.0;
    let hedeg: f64 = 20.0;
    let xnp: f64 = 3.0;
    let signoise: f64 = 0.001;
    let ts: f64 = 0.1;
    let tau: f64 = 0.5;
    let apn: i32 = 2;
    let tf: f64 = 10.0;

    let ts2 = ts * ts;
    let ts3 = ts2 * ts;
    let ts4 = ts3 * ts;
    let ts5 = ts4 * ts;
    let phin = xnt * xnt / tf;

    let mut array_tgos = Vec::new();
    let mut array_y = Vec::new();

    let mut tstart = 0.1;
    while tstart <= 10.0 {
        let mut y = yic;
        let mut yd: f64 = 0.0;
        let mut t: f64 = 0.0;
        let h: f64 = 0.01;
        let mut s: f64 = 0.0;

        let rtm = vc * tf;
        let sigpos = rtm * signoise;
        let sign2 = sigpos.powi(2);

        let mut p11 = sign2;
        let mut p12: f64 = 0.0;
        let mut p13: f64 = 0.0;
        let mut p22 = (vm * hedeg / 57.3).powi(2);
        let mut p23: f64 = 0.0;
        let mut p33 = xnt * xnt;

        let mut yh: f64 = 0.0;
        let mut ydh: f64 = 0.0;
        let mut xnth: f64 = 0.0;
        let mut xnc: f64 = 0.0;
        let mut xnl: f64 = 0.0;
        let mut xlam: f64;

        while t <= tf - 1e-5 {
            let yold = y;
            let ydold = yd;
            let xnlold = xnl;

            // First derivative evaluation
            let tgo = tf - t + 0.00001;
            let rtm = vc * tgo;
            let _xlam = y / (vc * tgo);
            let _xlamd = (rtm * yd + y * vc) / (rtm * rtm);
            let xnld = (xnc - xnl) / tau;
            let ydd = if t > tstart { xnt - xnl } else { 0.0 };

            // Euler step
            y += h * yd;
            yd += h * ydd;
            xnl += h * xnld;
            t += h;

            // Second derivative for RK2
            let tgo = tf - t + 0.00001;
            let rtm = vc * tgo;
            xlam = y / (vc * tgo); // Save XLAM for use in sampling block (matches MATLAB)
            let _xlamd = (rtm * yd + y * vc) / (rtm * rtm);
            let xnld = (xnc - xnl) / tau;
            let ydd = if t > tstart { xnt - xnl } else { 0.0 };

            // RK2 averaging
            y = 0.5 * (yold + y + h * yd);
            yd = 0.5 * (ydold + yd + h * ydd);
            xnl = 0.5 * (xnlold + xnl + h * xnld);

            s += h;

            if s >= ts - 1e-5 {
                s = 0.0;
                let tgo = tf - t + 0.000001;
                let rtm = vc * tgo;

                // Apply blind range - large measurement noise when TGO < TBLIND
                let sigpos = if tgo >= tblind {
                    rtm * signoise
                } else {
                    9999999999.0
                };
                let sign2 = sigpos.powi(2);

                // Kalman filter propagation
                let mut m11 = p11 + ts * p12 + 0.5 * ts2 * p13 + ts * (p12 + ts * p22 + 0.5 * ts2 * p23);
                m11 = m11 + 0.5 * ts2 * (p13 + ts * p23 + 0.5 * ts2 * p33) + ts5 * phin / 20.0;
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

                let xkblind = if tgo >= tblind { 1.0 } else { 0.0 };
                // Use XLAM from RK2 loop (not recomputed) to match MATLAB
                let ystar = rtm * xlam * xkblind;
                let res = ystar - yh - ts * ydh - 0.5 * ts * ts * (xnth - xnl);
                yh = k1 * res + yh + ts * ydh + 0.5 * ts * ts * (xnth - xnl);
                ydh = k2 * res + ydh + ts * (xnth - xnl);
                xnth += k3 * res;
                let xlamdh = (yh + ydh * tgo) / (vc * tgo * tgo);

                if apn == 0 {
                    let c1 = xnp / tgo.powi(2);
                    let c2 = xnp / tgo;
                    xnc = c1 * yh + c2 * ydh;
                } else if apn == 1 {
                    xnc = xnp * vc * xlamdh + 0.5 * xnp * xnth;
                } else {
                    let x = tgo / tau;
                    let top = 6.0 * x * x * ((-x).exp() - 1.0 + x);
                    let bot1 = 2.0 * x.powi(3) + 3.0 + 6.0 * x - 6.0 * x * x;
                    let bot2 = -12.0 * x * (-x).exp() - 3.0 * (-2.0 * x).exp();
                    let xnpp = top / (0.0001 + bot1 + bot2);
                    let xnew = xnpp * xnl * ((-x).exp() + x - 1.0) / (x * x);
                    xnc = xnpp * vc * xlamdh + 0.5 * xnpp * xnth - xnew;
                }

                if qzero == 1 && tgo < tblind {
                    xnc = 0.0;
                }
            }
        }

        let tgos = tf - tstart;
        array_tgos.push(tgos);
        array_y.push(y);

        tstart += 0.1;
    }

    Results {
        tgos: array_tgos,
        y: array_y,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c27l5_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.tgos.clone(),
        results.y.clone(),
    ])?;

    let plot_file = format!("{}/c27l5_miss.png", output_dir);
    let config = PlotConfig::new("Miss for Various TGO with Blind Range")
        .with_labels("Time to go at which maneuver occurs (S)", "Miss (Ft)");

    let series = vec![
        Series::new(results.tgos.clone(), results.y.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C27L5: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c27l5_runs() {
        let results = run();
        assert!(!results.tgos.is_empty());
    }
}
