//! Chapter 27, Lesson 6: Target Maneuver Miss with Pre-computed Gains and Blind Range
//!
//! Similar to L4 but includes blind range (TBLIND).

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub tp: Vec<f64>,
    pub xmnt: Vec<f64>,
    pub k1: Vec<f64>,
}

/// Run the C27L6 simulation
pub fn run() -> Results {
    let xnt: f64 = 96.6;
    let xnp: f64 = 3.0;
    let tf: f64 = 10.0;
    let ts: f64 = 0.1;
    let qzero: i32 = 0;
    let tblind: f64 = 0.5;
    let apn: i32 = 2;
    let vm: f64 = 3000.0;
    let vc: f64 = 4000.0;
    let hedeg: f64 = 20.0;
    let signoise: f64 = 0.001;

    let ts2 = ts * ts;
    let ts3 = ts2 * ts;
    let ts4 = ts3 * ts;
    let ts5 = ts4 * ts;
    let phin = xnt * xnt / tf;

    let rtm = vc * tf;
    let sigpos = rtm * signoise;
    let sign2 = sigpos.powi(2);

    let mut p11 = sign2;
    let mut p12: f64 = 0.0;
    let mut p13: f64 = 0.0;
    let mut p22 = (vm * hedeg / 57.3).powi(2);
    let mut p23: f64 = 0.0;
    let mut p33 = xnt * xnt;

    // Pre-compute Kalman gains with blind range
    // MATLAB initializes C=1 and does C=C+1 before storing, so U(1) is never set (implicitly 0)
    // We match this by starting with a zero
    let mut u_vec = vec![0.0];
    let mut v_vec = vec![0.0];
    let mut w_vec = vec![0.0];

    let mut t = ts;
    while t <= tf {
        let tgo = tf - t + 0.000001;
        let rtm = vc * tgo;

        // Apply blind range
        let sigpos = if tgo >= tblind {
            rtm * signoise
        } else {
            9999999999.0
        };
        let sign2 = sigpos.powi(2);

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

        // Store gains twice per MATLAB code (lines 55-61)
        u_vec.push(k1);
        v_vec.push(k2);
        w_vec.push(k3);
        u_vec.push(k1);
        v_vec.push(k2);
        w_vec.push(k3);

        t += ts;
    }

    // Reverse gains for time-to-go indexing
    let icount = (tf / ts).round() as usize;
    let mut e_vec = vec![0.0; icount];
    let mut f_vec = vec![0.0; icount];
    let mut g_vec = vec![0.0; icount];

    for i in 0..(icount - 1) {
        let rev = icount - 1 - i - 1;
        if i < u_vec.len() {
            e_vec[rev] = u_vec[i];
            f_vec[rev] = v_vec[i];
            g_vec[rev] = w_vec[i];
        }
    }

    // Main simulation
    let tap: f64 = 0.5;
    let h: f64 = 0.01;

    let mut tp = 0.00001;
    let mut s: f64 = 0.0;

    let mut x1: f64 = 0.0;
    let mut x2: f64 = 0.0;
    let mut x3: f64 = 1.0;
    let mut x4: f64 = 0.0;
    let mut x5: f64 = 0.0;
    let mut x6: f64 = 0.0;

    let mut y1old: f64 = 0.0;
    let mut y2old: f64 = 0.0;
    let mut y3old: f64 = 0.0;
    let mut y4old: f64 = 0.0;
    let mut y6old: f64 = 0.0;
    let mut y7old: f64 = 0.0;

    let mut y1new: f64 = 0.0;
    let mut y7new: f64 = 0.0;

    let mut array_tp = Vec::new();
    let mut array_xmnt = Vec::new();
    let mut array_k1 = Vec::new();

    let mut idx: usize = 0;

    while tp <= tf - 1e-5 {
        let x1old = x1;
        let x2old = x2;
        let x3old = x3;
        let x4old = x4;
        let x5old = x5;
        let x6old = x6;

        // First derivative evaluation
        let tgo = tp;
        let x1d = x2;
        let x2d = x3 + y1new / (vc * tgo);
        let x3d = y1new / (vc * tgo * tgo);
        let x4d = (x5 + y7new + x6) / tap;
        let x5d = -x4d;
        let x6d = -x2;

        // Euler step
        x1 += h * x1d;
        x2 += h * x2d;
        x3 += h * x3d;
        x4 += h * x4d;
        x5 += h * x5d;
        x6 += h * x6d;
        tp += h;

        // Second derivative for RK2
        let tgo = tp;
        let x1d = x2;
        let x2d = x3 + y1new / (vc * tgo);
        let x3d = y1new / (vc * tgo * tgo);
        let x4d = (x5 + y7new + x6) / tap;
        let x5d = -x4d;
        let x6d = -x2;

        // RK2 averaging
        x1 = (x1old + x1) / 2.0 + 0.5 * h * x1d;
        x2 = (x2old + x2) / 2.0 + 0.5 * h * x2d;
        x3 = (x3old + x3) / 2.0 + 0.5 * h * x3d;
        x4 = (x4old + x4) / 2.0 + 0.5 * h * x4d;
        x5 = (x5old + x5) / 2.0 + 0.5 * h * x5d;
        x6 = (x6old + x6) / 2.0 + 0.5 * h * x6d;

        s += h;

        if s > ts - 0.0001 {
            s = 0.0;

            let k1 = if idx < e_vec.len() { e_vec[idx] } else { 0.0 };
            let k2 = if idx < f_vec.len() { f_vec[idx] } else { 0.0 };
            let k3 = if idx < g_vec.len() { g_vec[idx] } else { 0.0 };
            idx += 1;

            let (c1, c2, c3, c4) = if apn == 0 {
                (xnp / tp.powi(2), xnp / tp, 0.0, 0.0)
            } else if apn == 1 {
                (xnp / tp.powi(2), xnp / tp, 0.5 * xnp, 0.0)
            } else {
                let x = tp / tap;
                let top = 6.0 * x * x * ((-x).exp() - 1.0 + x);
                let bot1 = 2.0 * x.powi(3) + 3.0 + 6.0 * x - 6.0 * x * x;
                let bot2 = -12.0 * x * (-x).exp() - 3.0 * (-2.0 * x).exp();
                let xnpp = top / (0.0001 + bot1 + bot2);
                (
                    xnpp / tp.powi(2),
                    xnpp / tp,
                    0.5 * xnpp,
                    -xnpp * ((-x).exp() + x - 1.0) / (x * x),
                )
            };

            let _tgo = tp;

            // XKBLIND for input (QZERO check)
            let xkblind = if tp < tblind && qzero == 1 { 0.0 } else { 1.0 };
            let temp1 = xkblind * (x4 - y6old);
            let temp2 = c1 * temp1 + y2old;
            let temp3 = c2 * temp1 + y3old;
            let temp4 = c3 * temp1 + y4old;
            let temp5 = k1 * temp2 + k2 * temp3 + k3 * temp4;

            // XKBLIND1 for output
            let xkblind1 = if tp >= tblind { 1.0 } else { 0.0 };
            y1new = y1old + temp5 * vc * tp * xkblind1;
            let y2new = temp2 - temp5;
            let y3new = temp3 + ts * y2new;
            let y4new = temp4 + ts * temp3 + 0.5 * ts * ts * y2new;
            let y5 = -(ts * temp3 + 0.5 * ts * ts * y2new);
            y7new = y7old + c4 * temp1 + y5;
            let y6new = x4;

            let xmnt = xnt * x1;

            y1old = y1new;
            y2old = y2new;
            y3old = y3new;
            y4old = y4new;
            y6old = y6new;
            y7old = y7new;

            array_tp.push(tp);
            array_xmnt.push(xmnt);
            array_k1.push(k1);
        }
    }

    Results {
        tp: array_tp,
        xmnt: array_xmnt,
        k1: array_k1,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c27l6_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.tp.clone(),
        results.xmnt.clone(),
        results.k1.clone(),
    ])?;

    let plot_file = format!("{}/c27l6_miss.png", output_dir);
    let config = PlotConfig::new("Target Maneuver Miss with Blind Range")
        .with_labels("Time to go at which maneuver occurs (S)", "Miss (Ft)");

    let series = vec![
        Series::new(results.tp.clone(), results.xmnt.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C27L6: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c27l6_runs() {
        let results = run();
        assert!(!results.tp.is_empty());
    }
}
