//! Chapter 30, Lesson 3: Kalman Filter with Frequency Estimation
//!
//! Implements a Kalman filter that estimates target weave frequency
//! in addition to acceleration states.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;
use rand::SeedableRng;
use rand_distr::{Distribution, Normal};

pub struct Results {
    pub time: Vec<f64>,
    pub ytddg: Vec<f64>,
    pub ytddhg: Vec<f64>,
    pub ytdddg: Vec<f64>,
    pub ytdddhg: Vec<f64>,
    pub w: Vec<f64>,
    pub wh: Vec<f64>,
}

/// PROJECT function - propagates state estimates forward
#[allow(clippy::too_many_arguments)]
fn project(ts: f64, yph: f64, ydph: f64, ytddph: f64, ytdddph: f64, hp: f64, xnlp: f64, wph: f64) -> (f64, f64, f64, f64) {
    let mut t: f64 = 0.0;
    let mut y = yph;
    let mut yd = ydph;
    let mut ytdd = ytddph;
    let mut ytddd = ytdddph;
    let w = wph;
    let xnl = xnlp;
    let h = hp;

    while t <= ts - 0.0001 {
        let ytdddd = -w * w * ytdd;
        ytddd += h * ytdddd;
        ytdd += h * ytddd;
        let ydd = ytdd - xnl;
        yd += h * ydd;
        y += h * yd;
        t += h;
    }

    (y, yd, ytdd, ytddd)
}

/// Run the C30L3 simulation
pub fn run() -> Results {
    run_with_seed(None)
}

/// Run with optional seed for reproducibility
pub fn run_with_seed(seed: Option<u64>) -> Results {
    let phis2: f64 = 0.0;
    let xnt: f64 = 96.6;
    let w: f64 = 2.0;
    let _phasedeg: f64 = 0.0;
    let sigrin: f64 = 0.0001;
    let siggl: f64 = 0.0;
    let srn: f64 = 0.0;
    let ra: f64 = 21000.0;
    let whic: f64 = -1.0;
    let ts: f64 = 0.01;
    let tf: f64 = 10.0;
    let phis1 = w * w * xnt * xnt / tf;
    let qperfect: bool = false;
    let vc: f64 = 9000.0;
    let xnp: f64 = 3.0;
    let xnclim: f64 = 9999999.0;
    let apn: i32 = 4;
    let tau: f64 = 0.5;
    let hedeg: f64 = 0.0;
    let vm: f64 = 3000.0;
    let order: usize = 5;

    let mut rng: rand::rngs::StdRng = match seed {
        Some(s) => rand::rngs::StdRng::seed_from_u64(s),
        None => rand::rngs::StdRng::from_entropy(),
    };
    let normal = Normal::new(0.0, 1.0).unwrap();

    let mut tgo: f64;
    let mut t: f64 = 0.0;
    let mut s: f64 = 0.0;
    let mut y: f64 = 0.0;
    let mut yd = -xnt / w - vm * hedeg / 57.3;
    let mut ytdd = xnt * (w * t).sin();
    let mut ytddd = xnt * w * (w * t).cos();
    let mut xnc: f64 = 0.0;
    let mut xnl: f64 = 0.0;
    let h: f64 = 0.001;
    let hp: f64 = 0.001;
    let ts2 = ts * ts;
    let _ts3 = ts2 * ts;

    let mut wh = whic;
    let (mut yh, mut ydh, mut ytddh, mut ytdddh) = if qperfect {
        (y, yd, ytdd, ytddd)
    } else {
        (0.0, 0.0, 0.0, 0.0)
    };
    if qperfect {
        wh = w;
    }

    // Initialize P matrix (5x5)
    let mut p = [[0.0f64; 5]; 5];
    let mut rtm = vc * tf;
    let signoise = (sigrin * sigrin + (siggl / rtm).powi(2) + (srn * rtm * rtm / (ra * ra)).powi(2)).sqrt();
    let ynoise = signoise * rtm;
    p[0][0] = ynoise * ynoise;
    p[1][1] = (vm * 20.0 / 57.3).powi(2);
    p[2][2] = xnt * xnt;
    p[3][3] = (w * xnt).powi(2);
    p[4][4] = w * w;

    let hmat = [1.0, 0.0, 0.0, 0.0, 0.0];

    let mut array_t = Vec::new();
    let mut array_ytddg = Vec::new();
    let mut array_ytddhg = Vec::new();
    let mut array_ytdddg = Vec::new();
    let mut array_ytdddhg = Vec::new();
    let mut array_w = Vec::new();
    let mut array_wh = Vec::new();

    while t <= tf - 0.0001 {
        let yold = y;
        let ydold = yd;
        let xnlold = xnl;

        // First derivative evaluation
        ytdd = xnt * (w * t).sin();
        let _tgo1 = tf - t + 0.00001;
        let xnld = (xnc - xnl) / tau;
        let ydd = ytdd - xnl;

        // Euler step
        y += h * yd;
        yd += h * ydd;
        xnl += h * xnld;
        t += h;

        // Second derivative evaluation
        ytdd = xnt * (w * t).sin();
        tgo = tf - t + 0.00001;
        let xnld = (xnc - xnl) / tau;
        let ydd = ytdd - xnl;

        // RK2 averaging
        y = 0.5 * (yold + y + h * yd);
        yd = 0.5 * (ydold + yd + h * ydd);
        xnl = 0.5 * (xnlold + xnl + h * xnld);

        s += h;
        if s >= ts - 0.00001 {
            s = 0.0;
            ytdd = xnt * (w * t).sin();
            ytddd = xnt * w * (w * t).cos();

            // Build PHI matrix
            let mut phi = [[0.0f64; 5]; 5];
            phi[0][0] = 1.0;
            phi[0][1] = ts;
            phi[1][1] = 1.0;
            phi[1][2] = ts;
            phi[2][2] = 1.0;
            phi[2][3] = ts;
            phi[3][2] = -wh * wh * ts;
            phi[3][3] = 1.0;
            phi[3][4] = -2.0 * wh * ytddh * ts;
            phi[4][4] = 1.0;

            // Build Q matrix
            let mut q = [[0.0f64; 5]; 5];
            q[2][2] = phis1 * ts * ts * ts / 3.0;
            q[2][3] = phis1 * ts * ts / 2.0;
            q[3][2] = q[2][3];
            q[3][3] = 4.0 * wh * wh * ytddh * ytddh * phis2 * ts * ts * ts / 3.0 + phis1 * ts;
            q[3][4] = -wh * ytddh * ts * ts * phis2;
            q[4][3] = q[3][4];
            q[4][4] = phis2 * ts;

            // Propagate covariance: M = PHI * P * PHI' + Q
            let mut phip = [[0.0f64; 5]; 5];
            for i in 0..order {
                for j in 0..order {
                    for k in 0..order {
                        phip[i][j] += phi[i][k] * p[k][j];
                    }
                }
            }
            let mut m = [[0.0f64; 5]; 5];
            for i in 0..order {
                for j in 0..order {
                    for k in 0..order {
                        m[i][j] += phip[i][k] * phi[j][k];
                    }
                    m[i][j] += q[i][j];
                }
            }

            // Kalman gain calculation
            rtm = vc * tgo;
            let signoise = (sigrin * sigrin + (siggl / rtm).powi(2) + (srn * rtm * rtm / (ra * ra)).powi(2)).sqrt();
            let ynoise = signoise * rtm;
            let rmat = ynoise * ynoise;

            let mut mht = [0.0f64; 5];
            for i in 0..order {
                for j in 0..order {
                    mht[i] += m[i][j] * hmat[j];
                }
            }
            let mut hmht: f64 = 0.0;
            for j in 0..order {
                hmht += hmat[j] * mht[j];
            }
            let hmhtr = hmht + rmat;
            let hmhtrinv = 1.0 / hmhtr;
            let mut gain = [0.0f64; 5];
            for i in 0..order {
                gain[i] = mht[i] * hmhtrinv;
            }

            // Update covariance
            for i in 0..order {
                for j in 0..order {
                    let mut ikh_m = 0.0f64;
                    for k in 0..order {
                        let ikh_ik = if i == k { 1.0 } else { 0.0 } - gain[i] * hmat[k];
                        ikh_m += ikh_ik * m[k][j];
                    }
                    p[i][j] = ikh_m;
                }
            }

            rtm = vc * tgo;
            let xlam = y / rtm;
            let xnoise = signoise * normal.sample(&mut rng);
            let xlams = xlam + xnoise;

            let (yb, ydb, ytddb, ytdddb) = project(ts, yh, ydh, ytddh, ytdddh, hp, xnl, wh);
            let res = rtm * xlams - yb;

            yh = yb + gain[0] * res;
            ydh = ydb + gain[1] * res;
            ytddh = ytddb + gain[2] * res;
            ytdddh = ytdddb + gain[3] * res;
            wh += gain[4] * res;

            let mut xnc_new = if apn == 0 {
                xnp * (yh + ydh * tgo) / (tgo * tgo)
            } else if apn == 1 {
                xnp * (yh + ydh * tgo + 0.5 * ytddh * tgo * tgo) / (tgo * tgo)
            } else if apn == 2 {
                let xs = tgo / tau;
                let top = 6.0 * xs * xs * ((-xs).exp() - 1.0 + xs);
                let bot1 = 2.0 * xs.powi(3) + 3.0 + 6.0 * xs - 6.0 * xs * xs;
                let bot2 = -12.0 * xs * (-xs).exp() - 3.0 * (-2.0 * xs).exp();
                let xnpp = top / (0.0001 + bot1 + bot2);
                let c1 = xnpp / (tgo * tgo);
                let c2 = xnpp / tgo;
                let c3 = 0.5 * xnpp;
                let c4 = -xnpp * ((-xs).exp() + xs - 1.0) / (xs * xs);
                c1 * yh + c2 * ydh + c3 * ytddh + c4 * xnl
            } else if apn == 3 {
                let xp = wh * tgo;
                xnp * (yh + ydh * tgo) / (tgo * tgo) + xnp * ytddh * (1.0 - xp.cos()) / (xp * xp)
                    + xnp * ytdddh * (xp - xp.sin()) / (xp * xp * wh)
            } else {
                let xs = tgo / tau;
                let top = 6.0 * xs * xs * ((-xs).exp() - 1.0 + xs);
                let bot1 = 2.0 * xs.powi(3) + 3.0 + 6.0 * xs - 6.0 * xs * xs;
                let bot2 = -12.0 * xs * (-xs).exp() - 3.0 * (-2.0 * xs).exp();
                let xnpp = top / (0.0001 + bot1 + bot2);
                let c1 = xnpp / (tgo * tgo);
                let c2 = xnpp / tgo;
                let c3 = xnpp * (1.0 - (wh * tgo).cos()) / (wh * wh * tgo * tgo);
                let c4 = -xnpp * ((-xs).exp() + xs - 1.0) / (xs * xs);
                let c5 = xnpp * (wh * tgo - (wh * tgo).sin()) / (wh * wh * wh * tgo * tgo);
                c1 * yh + c2 * ydh + c3 * ytddh + c4 * xnl + c5 * ytdddh
            };

            if xnc_new > xnclim {
                xnc_new = xnclim;
            }
            if xnc_new < -xnclim {
                xnc_new = -xnclim;
            }
            xnc = xnc_new;

            let ytddg = ytdd / 32.2;
            let ytddhg = ytddh / 32.2;
            let ytdddg = ytddd / 32.2;
            let ytdddhg = ytdddh / 32.2;

            array_t.push(t);
            array_ytddg.push(ytddg);
            array_ytddhg.push(ytddhg);
            array_ytdddg.push(ytdddg);
            array_ytdddhg.push(ytdddhg);
            array_w.push(w);
            array_wh.push(wh);
        }
    }

    Results {
        time: array_t,
        ytddg: array_ytddg,
        ytddhg: array_ytddhg,
        ytdddg: array_ytdddg,
        ytdddhg: array_ytdddhg,
        w: array_w,
        wh: array_wh,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c30l3_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.ytddg.clone(),
        results.ytddhg.clone(),
        results.ytdddg.clone(),
        results.ytdddhg.clone(),
        results.w.clone(),
        results.wh.clone(),
    ])?;

    let plot_file = format!("{}/c30l3_accel.png", output_dir);
    let config = PlotConfig::new("Acceleration and Estimate")
        .with_labels("Time (Sec)", "Acceleration (G)");

    let series = vec![
        Series::new(results.time.clone(), results.ytddg.clone())
            .with_color(plotters::prelude::BLUE)
            .with_label("Actual"),
        Series::new(results.time.clone(), results.ytddhg.clone())
            .with_color(plotters::prelude::RED)
            .with_label("Estimated"),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C30L3: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c30l3_runs() {
        let results = run_with_seed(Some(12345));
        assert!(!results.time.is_empty());
    }
}
