//! Chapter 9, Lesson 5: Optimal Guidance with Fifth-Order Binomial Filter
//!
//! Simulation with optimal proportional navigation and fifth-order binomial filter.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub tf: Vec<f64>,
    pub y: Vec<f64>,
}

/// Run the C9L5 simulation
pub fn run() -> Results {
    let vc: f64 = 4000.0;
    let xnt: f64 = 32.2;
    let yic: f64 = 0.0;
    let vm: f64 = 3000.0;
    let hedeg: f64 = 20.0;
    let _xnp: f64 = 3.0;
    let signoise: f64 = 0.001;
    let ts: f64 = 0.1;
    let tau: f64 = 0.5;
    let xlim: f64 = 644.0;
    let _r: f64 = 0.0;
    let ta: f64 = 5.0;
    let h: f64 = 0.01;

    let ts2 = ts * ts;
    let ts3 = ts2 * ts;
    let ts4 = ts3 * ts;
    let ts5 = ts4 * ts;

    let mut array_tf = Vec::new();
    let mut array_y = Vec::new();

    let mut tf_val: f64 = 0.1;
    while tf_val <= 10.0 {
        let mut y = yic;
        let mut yd = 0.0;
        let phin = xnt * xnt / tf_val;
        let rtm_init = vc * tf_val;
        let sigpos = rtm_init * signoise;
        let sign2_init = sigpos * sigpos;

        let mut p11 = sign2_init;
        let mut p12 = 0.0;
        let mut p13 = 0.0;
        let mut p22 = (vm * hedeg / 57.3).powi(2);
        let mut p23 = 0.0;
        let mut p33 = xnt * xnt;

        let mut t: f64 = 0.0;
        let mut s: f64 = 0.0;
        let mut yh = 0.0;
        let mut ydh = 0.0;
        let mut xnth = 0.0;
        let mut xnc = 0.0;
        let mut th = 0.0;
        let mut xnl = 0.0;
        let mut d = 0.0;
        let mut x4 = 0.0;
        let mut x5 = 0.0;
        let mut thh = 0.0;

        while t <= (tf_val - 1e-5) {
            let y_old = y;
            let yd_old = yd;
            let xnl_old = xnl;
            let th_old = th;
            let x4_old = x4;
            let x5_old = x5;
            let d_old = d;
            let thh_old = thh;

            // First derivative
            let tgo = tf_val - t + 0.00001;
            let _rtm = vc * tgo;
            let xlam = y / (vc * tgo);
            let eps = xlam - th - thh;
            let dd = 5.0 * eps / tau;
            let _xlams = eps + d;
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
            th += h * thd;
            x4 += h * x4d;
            x5 += h * x5d;
            d += h * dd;
            thh += h * thhd;
            t += h;

            // Second derivative
            let tgo = tf_val - t + 0.00001;
            let _rtm = vc * tgo;
            let xlam = y / (vc * tgo);
            let eps = xlam - th - thh;
            let dd = 5.0 * eps / tau;
            let xlams = eps + d;
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
            th = 0.5 * (th_old + th + h * thd);
            x4 = 0.5 * (x4_old + x4 + h * x4d);
            x5 = 0.5 * (x5_old + x5 + h * x5d);
            d = 0.5 * (d_old + d + h * dd);
            thh = 0.5 * (thh_old + thh + h * thhd);

            s += h;
            if s >= (ts - 1e-5) {
                s = 0.0;
                let tgo = tf_val - t + 0.000001;
                let rtm = vc * tgo;
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

                let xlamnoise = 0.0;
                let ystar = rtm * (xlams + xlamnoise);
                let res = ystar - yh - ts * ydh - 0.5 * ts * ts * (xnth - xnl);

                yh = k1 * res + yh + ts * ydh + 0.5 * ts * ts * (xnth - xnl);
                ydh = k2 * res + ydh + ts * (xnth - xnl);
                xnth += k3 * res;

                let xlamdh = (yh + ydh * tgo) / (vc * tgo * tgo);

                // Optimal guidance
                let x = tgo / tau;
                let top = 6.0 * x * x * ((-x).exp() - 1.0 + x);
                let bot1 = 2.0 * x * x * x + 3.0 + 6.0 * x - 6.0 * x * x;
                let bot2 = -12.0 * x * (-x).exp() - 3.0 * (-2.0 * x).exp();
                let xnpp = top / (0.0001 + bot1 + bot2);
                let xnew = xnpp * xnl * ((-x).exp() + x - 1.0) / (x * x);
                xnc = xnpp * vc * xlamdh + 0.5 * xnpp * xnth - xnew;

                if xnc > xlim {
                    xnc = xlim;
                } else if xnc < -xlim {
                    xnc = -xlim;
                }
            }
        }

        array_tf.push(tf_val);
        array_y.push(y);

        tf_val += 0.1;
    }

    Results {
        tf: array_tf,
        y: array_y,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c9l5_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.tf.clone(),
        results.y.clone(),
    ])?;

    let plot_file = format!("{}/c9l5_plot.png", output_dir);
    let config = PlotConfig::new("Miss for various flight times")
        .with_labels("Flight Time (S)", "Miss (Ft)");
    let series = vec![
        Series::new(results.tf.clone(), results.y.clone())
            .with_color(plotters::prelude::BLUE),
    ];
    line_plot(&plot_file, &config, &series).ok();

    println!("C9L5: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c9l5_runs() {
        let results = run();
        assert!(!results.tf.is_empty());
    }
}
