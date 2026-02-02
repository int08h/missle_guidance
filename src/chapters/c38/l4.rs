//! Chapter 38, Lesson 4: Kalman Filter with Trapezoidal Weave Target
//!
//! Single-run simulation using Kalman filter to estimate target acceleration
//! with trapezoidal weave target maneuver.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;
use rand::SeedableRng;
use rand_distr::{Distribution, Normal};

// Use MATLAB's PI value (3.1416) instead of std::f64::consts::PI
// to ensure numerical compatibility with the original code.
const PI: f64 = 3.1416;

pub struct Results {
    pub time: Vec<f64>,
    pub ytdd: Vec<f64>,
    pub xnth: Vec<f64>,
}

/// PROJECT6S function - project state forward for 6-state Kalman filter
#[allow(clippy::too_many_arguments)]
fn project6s(
    _tp: f64,
    ts: f64,
    yp: f64,
    ydp: f64,
    xnlp: f64,
    x1hp: f64,
    x2hp: f64,
    x3hp: f64,
    x4hp: f64,
    w: f64,
    hp: f64,
) -> (f64, f64, f64, f64, f64, f64) {
    let mut t: f64 = 0.0;
    let mut y = yp;
    let mut yd = ydp;
    let mut x1 = x1hp;
    let mut x2 = x2hp;
    let mut x3 = x3hp;
    let mut x4 = x4hp;
    let xnl = xnlp;
    let h = hp;

    while t <= ts - 0.0001 {
        let xnt = x1 + x3;
        let ydd = xnt - xnl;
        yd += h * ydd;
        y += h * yd;
        let x1d = x2;
        x1 += h * x1d;
        let x2d = -w * w * x1;
        x2 += h * x2d;
        let x3d = x4;
        x3 += h * x3d;
        let x4d = -9.0 * w * w * x3;
        x4 += h * x4d;
        t += h;
    }

    (y, yd, x1, x2, x3, x4)
}

/// Compute trapezoidal weave acceleration
fn trapezoidal_weave(t: f64, tstart: f64, xnt: f64, pz: f64, tr: f64, x: f64) -> f64 {
    if t < tstart {
        return 0.0;
    }

    let t_rel = t - tstart;
    let period_num = (t_rel / pz).floor() as i32;
    let tstar = t_rel - (period_num as f64) * pz;

    if tstar < tr / 2.0 {
        2.0 * xnt * tstar / tr
    } else if tstar < tr / 2.0 + x {
        xnt
    } else if tstar < 3.0 * tr / 2.0 + x {
        -2.0 * xnt * tstar / tr + 2.0 * xnt + 2.0 * xnt * x / tr
    } else if tstar < 3.0 * tr / 2.0 + 2.0 * x {
        -xnt
    } else {
        2.0 * xnt * tstar / tr - 4.0 * xnt - 4.0 * xnt * x / tr
    }
}

/// Run the C38L4 simulation
pub fn run() -> Results {
    run_with_seed(None)
}

/// Run with optional seed for reproducibility
#[allow(clippy::needless_range_loop)]
pub fn run_with_seed(seed: Option<u64>) -> Results {
    let _num_runs: usize = 100;
    let tau: f64 = 0.05;
    let order: usize = 6;
    let tstart: f64 = 0.0;
    let pz: f64 = 3.0;
    let tr: f64 = 0.5;
    let xl = pz / 2.0;
    let w = 2.0 * PI / pz;
    let _wh = w;
    let _wreal = w;
    let vc: f64 = 2.5 * 3280.0;
    let xnt: f64 = 322.0;
    let _xntreal = xnt;
    let ts: f64 = 0.01;
    let vm: f64 = 2.0 * 3280.0;
    let hedegfil: f64 = 20.0;
    let xnp: f64 = 3.0;
    let sigrin: f64 = 0.0001;
    let siggl: f64 = 0.0;
    let ra: f64 = 21000.0;
    let srn: f64 = 0.0;
    let tf: f64 = 10.0;
    let _qperfect: i32 = 0;
    let _phase = 0.0 / 57.3;
    let xlim: f64 = 1288.0;

    let _xntic = xnt;
    let x = pz / 2.0 - tr;
    let alf = PI * tr / (2.0 * xl);
    let phis = xnt * xnt / 6.0;

    let mut rng: rand::rngs::StdRng = match seed {
        Some(s) => rand::rngs::StdRng::seed_from_u64(s),
        None => rand::rngs::StdRng::from_entropy(),
    };
    let normal = Normal::new(0.0, 1.0).unwrap();

    let rtm = vc * tf;
    let signoise_init = (sigrin.powi(2) + (siggl / rtm).powi(2) + (srn * rtm * rtm / (ra * ra)).powi(2)).sqrt();
    let sigpos = rtm * signoise_init;
    let sign2 = sigpos * sigpos;

    // Initialize matrices
    let mut phi = [[0.0; 6]; 6];
    let mut p = [[0.0; 6]; 6];
    let mut qc = [[0.0; 6]; 6];
    let mut q = [[0.0; 6]; 6];
    let mut f = [[0.0; 6]; 6];
    let mut idnp = [[0.0; 6]; 6];

    for i in 0..order {
        idnp[i][i] = 1.0;
    }

    f[0][1] = 1.0;
    f[1][2] = 1.0;
    f[1][4] = 1.0;
    f[2][3] = 1.0;
    f[3][2] = -w * w;
    f[4][5] = 1.0;
    f[5][4] = -9.0 * w * w;

    // PHI = IDNP + F*TS + 0.5*F*F*TS*TS
    let mut ff = [[0.0; 6]; 6];
    for i in 0..order {
        for j in 0..order {
            for k in 0..order {
                ff[i][j] += f[i][k] * f[k][j];
            }
        }
    }

    for i in 0..order {
        for j in 0..order {
            phi[i][j] = idnp[i][j] + f[i][j] * ts + 0.5 * ff[i][j] * ts * ts;
        }
    }

    let tmpa = 4.0 * w * alf.sin() / (PI * alf);
    let tmpb = 4.0 * w * (3.0 * alf).sin() / (3.0 * PI * alf);
    qc[3][3] = phis * tmpa * tmpa;
    qc[3][5] = phis * tmpa * tmpb;
    qc[5][3] = phis * tmpa * tmpb;
    qc[5][5] = phis * tmpb * tmpb;

    // Q = PHI * QC * PHI' * TS
    let mut phiqc = [[0.0; 6]; 6];
    for i in 0..order {
        for j in 0..order {
            for k in 0..order {
                phiqc[i][j] += phi[i][k] * qc[k][j];
            }
        }
    }
    for i in 0..order {
        for j in 0..order {
            for k in 0..order {
                q[i][j] += phiqc[i][k] * phi[j][k];
            }
            q[i][j] *= ts;
        }
    }

    p[0][0] = sign2;
    p[1][1] = (vm * hedegfil / 57.3).powi(2);
    p[2][2] = 322.0_f64.powi(2);
    p[3][3] = (w * 322.0).powi(2);
    p[4][4] = 322.0_f64.powi(2);
    p[5][5] = (w * 322.0).powi(2);

    let hmat = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];

    let mut t: f64 = 0.0;
    let h: f64 = 0.001;
    let mut s: f64 = 0.0;
    let mut xnc: f64 = 0.0;
    let mut xnl: f64 = 0.0;
    let mut _ytdd: f64 = 0.0;
    let mut _ytddd: f64 = 0.0;
    let mut yh: f64 = 0.0;
    let mut ydh: f64 = 0.0;
    let mut x1h: f64 = 0.0;
    let mut x2h: f64 = 0.0;
    let mut x3h: f64 = 0.0;
    let mut x4h: f64 = 0.0;
    let mut y: f64 = 0.0;
    let mut yd: f64 = 0.0;

    let mut array_t = Vec::new();
    let mut array_ytdd = Vec::new();
    let mut array_xnth = Vec::new();

    while t <= tf - 0.0001 {
        let yold = y;
        let ydold = yd;
        let xnlold = xnl;

        // First derivative evaluation
        let tgo = tf - t + 0.000001;
        let _rtm = vc * tgo;
        let _xlam = y / (vc * tgo);
        let ytdd = trapezoidal_weave(t, tstart, xnt, pz, tr, x);
        let xnld = (xnc - xnl) / tau;
        let ydd = ytdd - xnl;

        // Euler step
        y += h * yd;
        yd += h * ydd;
        xnl += h * xnld;
        t += h;

        // Second derivative for RK2
        let _tgo = tf - t + 0.000001;
        let _rtm = vc * _tgo;
        let _xlam = y / (vc * _tgo);
        let ytdd = trapezoidal_weave(t, tstart, xnt, pz, tr, x);
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
            let rtm = vc * tgo;
            let signoise = (sigrin.powi(2) + (siggl / rtm).powi(2)
                + (srn * rtm * rtm / (ra * ra)).powi(2)).sqrt();
            let sigpos = rtm * signoise;
            let sign2 = sigpos * sigpos;
            let rmat = sign2;

            // M = PHI * P * PHI' + Q
            let mut phip = [[0.0; 6]; 6];
            for i in 0..order {
                for j in 0..order {
                    for k in 0..order {
                        phip[i][j] += phi[i][k] * p[k][j];
                    }
                }
            }
            let mut m = [[0.0; 6]; 6];
            for i in 0..order {
                for j in 0..order {
                    for k in 0..order {
                        m[i][j] += phip[i][k] * phi[j][k];
                    }
                    m[i][j] += q[i][j];
                }
            }

            // K = M * HMAT' / (HMAT * M * HMAT' + RMAT)
            let hmh = m[0][0] + rmat;
            let mut k = [0.0; 6];
            for i in 0..order {
                k[i] = m[i][0] / hmh;
            }

            // P = (IDNP - K*HMAT) * M
            let mut kh = [[0.0; 6]; 6];
            for i in 0..order {
                kh[i][0] = k[i] * hmat[0];
            }
            let mut idkh = [[0.0; 6]; 6];
            for i in 0..order {
                for j in 0..order {
                    idkh[i][j] = idnp[i][j] - kh[i][j];
                }
            }
            let mut new_p = [[0.0; 6]; 6];
            for i in 0..order {
                for j in 0..order {
                    for kk in 0..order {
                        new_p[i][j] += idkh[i][kk] * m[kk][j];
                    }
                }
            }
            p = new_p;

            let (yb, ydb, x1b, x2b, x3b, x4b) = project6s(
                t, ts, yh, ydh, xnl, x1h, x2h, x3h, x4h, w, h,
            );

            let xlam = y / (vc * tgo);
            let xlamnoise = signoise * normal.sample(&mut rng);
            let ystar = rtm * (xlam + xlamnoise);
            let res = ystar - yb;
            yh = yb + k[0] * res;
            ydh = ydb + k[1] * res;
            x1h = x1b + k[2] * res;
            x2h = x2b + k[3] * res;
            x3h = x3b + k[4] * res;
            x4h = x4b + k[5] * res;

            let xnth = x1h + x3h;
            xnc = xnp * (yh + ydh * tgo + 0.5 * xnth * tgo * tgo) / (tgo * tgo);
            if xnc > xlim {
                xnc = xlim;
            }
            if xnc < -xlim {
                xnc = -xlim;
            }

            let ytdd = trapezoidal_weave(t, tstart, xnt, pz, tr, x);
            array_t.push(t);
            array_ytdd.push(ytdd);
            array_xnth.push(xnth);
        }
    }

    Results {
        time: array_t,
        ytdd: array_ytdd,
        xnth: array_xnth,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c38l4_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.ytdd.clone(),
        results.xnth.clone(),
    ])?;

    let plot_file = format!("{}/c38l4_acc.png", output_dir);
    let config = PlotConfig::new("Target Acceleration Estimation")
        .with_labels("Flight Time (S)", "Acceleration (ft/s^2)");

    let series = vec![
        Series::new(results.time.clone(), results.ytdd.clone())
            .with_color(plotters::prelude::BLUE)
            .with_label("Actual"),
        Series::new(results.time.clone(), results.xnth.clone())
            .with_color(plotters::prelude::RED)
            .with_label("Estimate"),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C38L4: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c38l4_runs() {
        let results = run_with_seed(Some(12345));
        assert!(!results.time.is_empty());
    }
}
