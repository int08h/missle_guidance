//! Chapter 30, Lesson 2: Kalman Filter with Singer Model
//!
//! Implements a Kalman filter using the Singer target acceleration model
//! for estimating target acceleration and jerk.

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
}

/// Run the C30L2 simulation
pub fn run() -> Results {
    run_with_seed(None)
}

/// Run with optional seed for reproducibility
pub fn run_with_seed(seed: Option<u64>) -> Results {
    let tau: f64 = 0.5;
    let apn: i32 = 0;
    let order: usize = 4;
    let mvr: i32 = 1;
    let vc: f64 = 9000.0;
    let w: f64 = 2.0;
    let wreal: f64 = 2.0;
    let wh: f64 = w;
    let xnt: f64 = 96.6;
    let xntreal: f64 = 96.6;
    let ts: f64 = 0.01;
    let yic: f64 = 0.0;
    let vm: f64 = 3000.0;
    let hedeg: f64 = 0.0;
    let hedegfil: f64 = 20.0;
    let xnp: f64 = 3.0;
    let sigrin: f64 = 0.001;
    let siggl: f64 = 0.0;
    let ra: f64 = 21000.0;
    let srn: f64 = 0.0;
    let tf: f64 = 10.0;
    let qperfect: bool = false;

    let mut rng: rand::rngs::StdRng = match seed {
        Some(s) => rand::rngs::StdRng::seed_from_u64(s),
        None => rand::rngs::StdRng::from_entropy(),
    };
    let normal = Normal::new(0.0, 1.0).unwrap();

    let x_ts = wh * ts;
    let mut y = yic;
    let mut yd = -vm * hedeg / 57.3;
    let phis = wh * wh * xnt * xnt / tf;
    let rtm_init = vc * tf;
    let mut signoise = (sigrin * sigrin + (siggl / rtm_init).powi(2) + (srn * rtm_init * rtm_init / (ra * ra)).powi(2)).sqrt();
    let sigpos = rtm_init * signoise;
    let sign2 = sigpos * sigpos;

    // Initialize PHI matrix (4x4)
    let mut phi = [[0.0f64; 4]; 4];
    phi[0][0] = 1.0;
    phi[0][1] = ts;
    phi[0][2] = (1.0 - x_ts.cos()) / (wh * wh);
    phi[0][3] = (x_ts - x_ts.sin()) / (wh * wh * wh);
    phi[1][1] = 1.0;
    phi[1][2] = x_ts.sin() / wh;
    phi[1][3] = (1.0 - x_ts.cos()) / (wh * wh);
    phi[2][2] = x_ts.cos();
    phi[2][3] = x_ts.sin() / wh;
    phi[3][2] = -wh * x_ts.sin();
    phi[3][3] = x_ts.cos();

    // Initialize Q matrix
    let mut q = [[0.0f64; 4]; 4];
    q[0][0] = phis * (0.333 * x_ts.powi(3) - 2.0 * x_ts.sin() + 2.0 * x_ts * x_ts.cos() + 0.5 * x_ts - 0.25 * (2.0 * x_ts).sin()) / wh.powi(5);
    q[0][1] = phis * (0.5 * x_ts * x_ts - x_ts * x_ts.sin() + 0.5 * x_ts.sin() * x_ts.sin()) / wh.powi(4);
    q[1][0] = q[0][1];
    q[0][2] = phis * (x_ts.sin() - x_ts * x_ts.cos() - 0.5 * x_ts + 0.25 * (2.0 * x_ts).sin()) / wh.powi(3);
    q[2][0] = q[0][2];
    q[0][3] = phis * (x_ts.cos() + x_ts * x_ts.sin() - 0.5 * x_ts.sin() * x_ts.sin() - 1.0) / (wh * wh);
    q[3][0] = q[0][3];
    q[1][1] = phis * (1.5 * x_ts - 2.0 * x_ts.sin() + 0.25 * (2.0 * x_ts).sin()) / wh.powi(3);
    q[1][2] = phis * (1.0 - x_ts.cos() - 0.5 * x_ts.sin() * x_ts.sin()) / (wh * wh);
    q[2][1] = q[1][2];
    q[1][3] = phis * (x_ts.sin() - 0.5 * x_ts - 0.25 * (2.0 * x_ts).sin()) / wh;
    q[3][1] = q[1][3];
    q[2][2] = phis * (0.5 * x_ts - 0.25 * (2.0 * x_ts).sin()) / wh;
    q[2][3] = 0.5 * phis * x_ts.sin() * x_ts.sin();
    q[3][2] = q[2][3];
    q[3][3] = wh * phis * (0.5 * x_ts + 0.25 * (2.0 * x_ts).sin());

    // Initialize P matrix
    let mut p = [[0.0f64; 4]; 4];
    p[0][0] = sign2;
    p[1][1] = (vm * hedegfil / 57.3).powi(2);
    p[2][2] = xnt * xnt;
    p[3][3] = wh * wh * xnt * xnt;

    let hmat = [1.0, 0.0, 0.0, 0.0];

    let mut t: f64 = 0.0;
    let h: f64 = 0.001;
    let mut s: f64 = 0.0;
    let mut xnc: f64 = 0.0;
    let mut xnl: f64 = 0.0;
    let mut xlam: f64;

    let (mut yh, mut ydh, mut ytddh, mut ytdddh) = if qperfect {
        let ytdd = if mvr == 0 { xntreal } else { xntreal * (wreal * t).sin() };
        let ytddd = if mvr == 0 { 0.0 } else { xntreal * wreal * (wreal * t).cos() };
        (y, yd, ytdd, ytddd)
    } else {
        (0.0, 0.0, 0.0, 0.0)
    };

    let mut array_t = Vec::new();
    let mut array_ytddg = Vec::new();
    let mut array_ytddhg = Vec::new();
    let mut array_ytdddg = Vec::new();
    let mut array_ytdddhg = Vec::new();

    while t <= tf - 0.0001 {
        let yold = y;
        let ydold = yd;
        let xnlold = xnl;

        // First derivative evaluation
        let tgo = tf - t + 0.000001;
        let _rtm1 = vc * tgo;
        let _xlam1 = y / (vc * tgo);
        let ytdd = if mvr == 0 { xntreal } else { xntreal * (wreal * t).sin() };
        let xnld = (xnc - xnl) / tau;
        let ydd = ytdd - xnl;

        // Euler step
        y += h * yd;
        yd += h * ydd;
        xnl += h * xnld;
        t += h;

        // Second derivative evaluation
        let tgo = tf - t + 0.000001;
        let rtm: f64;
        xlam = y / (vc * tgo);
        let ytdd = if mvr == 0 { xntreal } else { xntreal * (wreal * t).sin() };
        let xnld = (xnc - xnl) / tau;
        let ydd = ytdd - xnl;

        // RK2 averaging
        y = 0.5 * (yold + y + h * yd);
        yd = 0.5 * (ydold + yd + h * ydd);
        xnl = 0.5 * (xnlold + xnl + h * xnld);

        s += h;
        if s >= ts - 0.00001 {
            s = 0.0;
            let tgo = tf - t + 0.000001;
            rtm = vc * tgo;
            signoise = (sigrin * sigrin + (siggl / rtm).powi(2) + (srn * rtm * rtm / (ra * ra)).powi(2)).sqrt();
            let sigpos = rtm * signoise;
            let sign2 = sigpos * sigpos;
            let rmat = sign2;

            // Propagate covariance: M = PHI * P * PHI' + Q
            let mut phip = [[0.0f64; 4]; 4];
            for i in 0..order {
                for j in 0..order {
                    for k in 0..order {
                        phip[i][j] += phi[i][k] * p[k][j];
                    }
                }
            }
            let mut m = [[0.0f64; 4]; 4];
            for i in 0..order {
                for j in 0..order {
                    for k in 0..order {
                        m[i][j] += phip[i][k] * phi[j][k];
                    }
                    m[i][j] += q[i][j];
                }
            }

            // Kalman gain: GAIN = M * H' / (H * M * H' + R)
            let mut mht = [0.0f64; 4];
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
            let mut gain = [0.0f64; 4];
            for i in 0..order {
                gain[i] = mht[i] * hmhtrinv;
            }

            // Update covariance: P = (I - K*H) * M
            let mut kh = [[0.0f64; 4]; 4];
            for i in 0..order {
                for j in 0..order {
                    kh[i][j] = gain[i] * hmat[j];
                }
            }
            for i in 0..order {
                for j in 0..order {
                    let _ikh = if i == j { 1.0 } else { 0.0 } - kh[i][j];
                    p[i][j] = 0.0;
                    for k in 0..order {
                        let ikh_k = if i == k { 1.0 } else { 0.0 } - kh[i][k];
                        p[i][j] += ikh_k * m[k][j];
                    }
                }
            }

            let (ytdd, ytddd) = if mvr == 0 {
                (xntreal, 0.0)
            } else {
                (xntreal * (wreal * t).sin(), xntreal * wreal * (wreal * t).cos())
            };

            let xlamnoise = signoise * normal.sample(&mut rng);
            let ystar = rtm * (xlam + xlamnoise);
            let res = ystar - yh - ts * ydh - (1.0 - x_ts.cos()) * ytddh / (wh * wh)
                - (x_ts - x_ts.sin()) * ytdddh / (wh * wh * wh) + 0.5 * ts * ts * xnl;

            yh = yh + ts * ydh + (1.0 - x_ts.cos()) * ytddh / (wh * wh)
                + (x_ts - x_ts.sin()) * ytdddh / (wh * wh * wh) + gain[0] * res - 0.5 * ts * ts * xnl;
            ydh = ydh + x_ts.sin() * ytddh / wh + (1.0 - x_ts.cos()) * ytdddh / (wh * wh)
                + gain[1] * res - ts * xnl;
            let ytddhnew = x_ts.cos() * ytddh + x_ts.sin() * ytdddh / wh + gain[2] * res;
            ytdddh = -wh * x_ts.sin() * ytddh + x_ts.cos() * ytdddh + gain[3] * res;
            ytddh = ytddhnew;

            if apn == 0 {
                xnc = xnp * (yh + ydh * tgo) / (tgo * tgo);
            } else if apn == 1 {
                xnc = xnp * (yh + ydh * tgo + 0.5 * ytddh * tgo * tgo) / (tgo * tgo);
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
                xnc = c1 * yh + c2 * ydh + c3 * ytddh + c4 * xnl;
            } else if apn == 3 {
                let xp = wh * tgo;
                xnc = xnp * (yh + ydh * tgo) / (tgo * tgo) + xnp * ytddh * (1.0 - xp.cos()) / (xp * xp)
                    + xnp * ytdddh * (xp - xp.sin()) / (xp * xp * wh);
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
                xnc = c1 * yh + c2 * ydh + c3 * ytddh + c4 * xnl + c5 * ytdddh;
            }

            let ytddg = ytdd / 32.2;
            let ytddhg = ytddh / 32.2;
            let ytdddg = ytddd / 32.2;
            let ytdddhg = ytdddh / 32.2;

            array_t.push(t);
            array_ytddg.push(ytddg);
            array_ytddhg.push(ytddhg);
            array_ytdddg.push(ytdddg);
            array_ytdddhg.push(ytdddhg);
        }
    }

    Results {
        time: array_t,
        ytddg: array_ytddg,
        ytddhg: array_ytddhg,
        ytdddg: array_ytdddg,
        ytdddhg: array_ytdddhg,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c30l2_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.ytddg.clone(),
        results.ytddhg.clone(),
        results.ytdddg.clone(),
        results.ytdddhg.clone(),
    ])?;

    let plot_file = format!("{}/c30l2_accel.png", output_dir);
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

    println!("C30L2: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c30l2_runs() {
        let results = run_with_seed(Some(12345));
        assert!(!results.time.is_empty());
    }
}
