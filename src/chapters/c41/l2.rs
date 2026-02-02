//! Chapter 41, Lesson 2: Monte Carlo Kalman Filter Miss Estimation
//!
//! Monte Carlo simulation with Kalman filtering for RMS miss distance.

use crate::save_data;
use rand::Rng;
use rand_distr::{Distribution, Normal};

pub struct Results {
    pub tf: Vec<f64>,
    pub rms: Vec<f64>,
    pub sp11: Vec<f64>,
    pub form: Vec<f64>,
}

/// Run the C41L2 simulation
pub fn run() -> Results {
    let vc: f64 = 5.0 * 3280.0;
    let xntic: f64 = 161.0;
    let vm: f64 = 3000.0;
    let hedeg: f64 = 20.0;
    let signoise: f64 = 10.0;
    let ts: f64 = 0.01;
    let tau: f64 = 0.2;
    let run_count: usize = 100;
    let amaxg: f64 = 99999999.0;
    let pz1: f64 = 0.0001;
    let amax = amaxg * 32.2;
    let phin = signoise * signoise * ts;

    let ts2 = ts * ts;
    let ts3 = ts2 * ts;
    let ts4 = ts3 * ts;
    let ts5 = ts4 * ts;

    let mut array_tf = Vec::new();
    let mut array_rms = Vec::new();
    let mut array_sp11 = Vec::new();
    let mut array_form = Vec::new();

    let mut rng = rand::thread_rng();
    let normal = Normal::new(0.0, 1.0).unwrap();

    let mut tf = 0.2;
    while tf <= 10.0 {
        let mut z1_sum: f64 = 0.0;
        let mut z_arr: Vec<f64> = Vec::with_capacity(run_count);

        let phis = xntic * xntic / tf;
        let mut sp11_last: f64 = 0.0;

        for _jj in 0..run_count {
            let sum_rand: f64 = rng.gen();
            let tstart = tf * sum_rand;
            let pz: f64 = rng.gen::<f64>() - 0.5;
            let coef: f64 = if pz > 0.0 { 1.0 } else { -1.0 };

            let mut y: f64 = 0.0;
            let mut yd: f64 = 0.0;
            let sigpos = signoise;
            let sign2 = sigpos * sigpos;
            let mut p11 = sign2;
            let mut p12: f64 = 0.0;
            let mut p13: f64 = 0.0;
            let mut p22 = (vm * hedeg / 57.3).powi(2);
            let mut p23: f64 = 0.0;
            let mut p33 = xntic * xntic;

            let mut t: f64 = 0.0;
            let h: f64 = 0.001;
            let mut s: f64 = 0.0;
            let mut yh: f64 = 0.0;
            let mut ydh: f64 = 0.0;
            let mut xnth: f64 = 0.0;
            let mut xnc: f64 = 0.0;
            let mut xnl: f64 = 0.0;

            while t <= tf - 1e-5 {
                let yold = y;
                let ydold = yd;
                let xnlold = xnl;

                // First derivative evaluation
                let xntc = if t < tstart { 0.0 } else { coef * xntic };
                let _tgo = tf - t + 0.00001;
                let xnld = (xnc - xnl) / tau;
                let xnt = xntc;
                let ydd = xnt - xnl;

                // Euler step
                y += h * yd;
                yd += h * ydd;
                xnl += h * xnld;
                t += h;

                // Second evaluation
                let xntc = if t < tstart { 0.0 } else { coef * xntic };
                let _tgo = tf - t + 0.00001;
                let xnld = (xnc - xnl) / tau;
                let xnt = xntc;
                let ydd = xnt - xnl;

                // RK2 averaging
                y = 0.5 * (yold + y + h * yd);
                yd = 0.5 * (ydold + yd + h * ydd);
                xnl = 0.5 * (xnlold + xnl + h * xnld);

                s += h;
                if s >= ts - 1e-5 {
                    s = 0.0;
                    let tgo = tf - t + 0.000001;
                    let _rtm = vc * tgo;
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

                    let k1 = m11 / (m11 + sign2);
                    let k2 = m12 / (m11 + sign2);
                    let k3 = m13 / (m11 + sign2);

                    p11 = (1.0 - k1) * m11;
                    p12 = (1.0 - k1) * m12;
                    p13 = (1.0 - k1) * m13;
                    p22 = -k2 * m12 + m22;
                    p23 = -k2 * m13 + m23;
                    p33 = -k3 * m13 + m33;

                    let ynoise = signoise * normal.sample(&mut rng);
                    let ystar = y + ynoise;
                    let res = ystar - yh - ts * ydh - 0.5 * ts * ts * (xnth - xnl);
                    yh = k1 * res + yh + ts * ydh + 0.5 * ts * ts * (xnth - xnl);
                    ydh = k2 * res + ydh + ts * (xnth - xnl);
                    xnth += k3 * res;

                    let x = tgo / tau;
                    let zem2h = yh + ydh * tgo - xnl * tau * tau * ((-x).exp() + x - 1.0)
                        + 0.5 * xnth * tgo * tgo;
                    let top = 6.0 * x * x * ((-x).exp() - 1.0 + x);
                    let bot1 = 2.0 * x * x * x + 3.0 + 6.0 * x - 6.0 * x * x;
                    let bot2 = -12.0 * x * (-x).exp() - 3.0 * (-2.0 * x).exp();
                    let xnpp = top / (pz1 + bot1 + bot2);
                    let _xnew = xnpp * xnl * ((-x).exp() + x - 1.0) / (x * x);
                    xnc = xnpp * zem2h / (tgo * tgo);
                    xnc = xnc.clamp(-amax, amax);
                }
            }

            sp11_last = p11.sqrt();
            z_arr.push(y);
            z1_sum += y;
        }

        // Compute statistics
        let _xmean = z1_sum / run_count as f64;
        let mut z2: f64 = 0.0;
        for z in &z_arr {
            z2 += z * z;
        }
        let rms = if run_count > 1 {
            (z2 / (run_count - 1) as f64).sqrt()
        } else {
            0.0
        };
        let form = (2.0 * phis.powf(0.16667) * phin.powf(0.8333)).sqrt();

        array_tf.push(tf);
        array_rms.push(rms);
        array_sp11.push(sp11_last);
        array_form.push(form);

        tf += 0.2;
    }

    Results {
        tf: array_tf,
        rms: array_rms,
        sp11: array_sp11,
        form: array_form,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c41l2_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.tf.clone(),
        results.rms.clone(),
        results.sp11.clone(),
        results.form.clone(),
    ])?;

    println!("C41L2: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c41l2_runs() {
        let results = run();
        assert!(!results.tf.is_empty());
    }
}
