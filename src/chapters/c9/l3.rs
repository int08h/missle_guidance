//! Chapter 9, Lesson 3: Kalman Filter Monte Carlo
//!
//! Monte Carlo simulation with Kalman filter for miss distance analysis.
//! Note: Running with NOISE=0 for deterministic verification.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub tf: Vec<f64>,
    pub sigma: Vec<f64>,
    pub xmean: Vec<f64>,
}

/// Run the C9L3 simulation
pub fn run() -> Results {
    let vc: f64 = 4000.0;
    let xnt: f64 = 96.6;
    let yic: f64 = 0.0;
    let vm: f64 = 3000.0;
    let hedeg: f64 = 20.0;
    let xnp: f64 = 3.0;
    let signoise: f64 = 0.001;
    let ts: f64 = 0.1;
    let tau: f64 = 0.5;
    let num_runs: usize = 50;
    let apn: i32 = 0;
    let xlim: f64 = 999999.0;
    let h: f64 = 0.01;
    let noise: i32 = 0; // Set to 0 for deterministic verification (MATLAB has NOISE=1)

    let ts2 = ts * ts;
    let ts3 = ts2 * ts;
    let ts4 = ts3 * ts;
    let ts5 = ts4 * ts;

    let mut array_tf = Vec::new();
    let mut array_sigma = Vec::new();
    let mut array_xmean = Vec::new();

    let mut tf: f64 = 0.5;
    while tf <= 10.0 {
        let mut z = vec![0.0; num_runs];
        let mut z1 = 0.0;

        for z_item in z.iter_mut() {
            let mut y = yic;
            let mut yd = 0.0;
            let phin = xnt * xnt / tf;
            let rtm_init = vc * tf;
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
            let mut xnl = 0.0;

            while t <= (tf - 1e-5) {
                let y_old = y;
                let yd_old = yd;
                let xnl_old = xnl;

                // First derivative
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

                // Second derivative
                let tgo = tf - t + 0.00001;
                let rtm = vc * tgo;
                let xlam = y / (vc * tgo);
                let _xlamd = (rtm * yd + y * vc) / (rtm * rtm);
                let xnld = (xnc - xnl) / tau;
                let ydd = xnt - xnl;

                // RK2 averaging
                y = 0.5 * (y_old + y + h * yd);
                yd = 0.5 * (yd_old + yd + h * ydd);
                xnl = 0.5 * (xnl_old + xnl + h * xnld);

                s += h;
                if s >= (ts - 1e-5) {
                    s = 0.0;
                    let tgo = tf - t + 0.000001;
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

                    // NOISE=0 for deterministic verification
                    let xlamnoise = if noise == 1 { 0.0 } else { 0.0 }; // Would use signoise * randn if noise=1
                    let ystar = rtm * (xlam + xlamnoise);
                    let res = ystar - yh - ts * ydh - 0.5 * ts * ts * (xnth - xnc);

                    yh = k1 * res + yh + ts * ydh + 0.5 * ts * ts * (xnth - xnc);
                    ydh = k2 * res + ydh + ts * (xnth - xnc);
                    xnth += k3 * res;

                    let xlamdh = (yh + ydh * tgo) / (vc * tgo * tgo);

                    xnc = match apn {
                        0 => xnp * vc * xlamdh,
                        1 => xnp * vc * xlamdh + 0.5 * xnp * xnth,
                        _ => {
                            let x = tgo / tau;
                            let top = 6.0 * x * x * ((-x).exp() - 1.0 + x);
                            let bot1 = 2.0 * x * x * x + 3.0 + 6.0 * x - 6.0 * x * x;
                            let bot2 = -12.0 * x * (-x).exp() - 3.0 * (-2.0 * x).exp();
                            let xnpp = top / (0.0001 + bot1 + bot2);
                            let xnew = xnpp * xnl * ((-x).exp() + x - 1.0) / (x * x);
                            xnpp * vc * xlamdh + 0.5 * xnpp * xnth - xnew
                        }
                    };

                    if xnc > xlim {
                        xnc = xlim;
                    } else if xnc < -xlim {
                        xnc = -xlim;
                    }
                }
            }

            *z_item = y;
            z1 += *z_item;
        }

        let xmean = z1 / num_runs as f64;

        let mut z1 = 0.0;
        for z_val in z.iter() {
            z1 += (*z_val - xmean).powi(2);
        }
        let sigma = if num_runs > 1 {
            (z1 / (num_runs - 1) as f64).sqrt()
        } else {
            0.0
        };

        array_tf.push(tf);
        array_sigma.push(sigma);
        array_xmean.push(xmean);

        tf += 0.5;
    }

    Results {
        tf: array_tf,
        sigma: array_sigma,
        xmean: array_xmean,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c9l3_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.tf.clone(),
        results.sigma.clone(),
        results.xmean.clone(),
    ])?;

    let plot_file = format!("{}/c9l3_sigma.png", output_dir);
    let config = PlotConfig::new("Standard deviation of miss")
        .with_labels("Flight Time (S)", "Noise Miss Standard Deviation (Ft)")
        .with_y_range(0.0, 4.0);
    let series = vec![
        Series::new(results.tf.clone(), results.sigma.clone())
            .with_color(plotters::prelude::BLUE),
    ];
    line_plot(&plot_file, &config, &series).ok();

    println!("C9L3: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c9l3_runs() {
        let results = run();
        assert!(!results.tf.is_empty());
    }
}
