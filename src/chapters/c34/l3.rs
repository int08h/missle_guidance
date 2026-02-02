//! Chapter 34, Lesson 3: Monte Carlo with Kalman Filter
//!
//! Monte Carlo simulation with random target maneuver timing
//! using Kalman filter estimation.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;
use rand::SeedableRng;
use rand_distr::{Distribution, Normal, Uniform};

pub struct Results {
    pub tf: Vec<f64>,
    pub rms: Vec<f64>,
}

/// PROJECT34 function - project state forward
#[allow(clippy::too_many_arguments)]
fn project34(
    _tp: f64,
    ts: f64,
    yhp: f64,
    ydhp: f64,
    xnlp: f64,
    xnthp: f64,
    xnu: f64,
    ipoisson: i32,
) -> (f64, f64, f64) {
    let mut t: f64 = 0.0;
    let mut yh = yhp;
    let mut ydh = ydhp;
    let mut xnt = xnthp;
    let xnl = xnlp;
    let h: f64 = 0.001;

    while t <= ts - 0.0001 {
        if ipoisson == 1 {
            let xntd = -2.0 * xnu * xnt;
            xnt += h * xntd;
            let yddh = xnt - xnl;
            ydh += h * yddh;
            yh += h * ydh;
        } else {
            let yddh = xnt - xnl;
            ydh += h * yddh;
            yh += h * ydh;
        }
        t += h;
    }

    (yh, ydh, xnt)
}

/// Run the C34L3 simulation
pub fn run() -> Results {
    run_with_seed(None)
}

/// Run with optional seed for reproducibility
#[allow(clippy::needless_range_loop)]
pub fn run_with_seed(seed: Option<u64>) -> Results {
    let num_runs: usize = 100;
    let order: usize = 3;
    let ipoisson: i32 = 1;
    let vc: f64 = 5.0 * 3280.0;
    let xntic: f64 = 161.0;
    let _yic: f64 = 0.0;
    let vm: f64 = 3000.0;
    let hedeg: f64 = 20.0;
    let _xnp: f64 = 3.0;
    let signoise: f64 = 0.0001;
    let ts: f64 = 0.01;
    let tau: f64 = 0.2;
    let amaxg: f64 = 10.0;
    let iconstant: i32 = 2;
    let _w: f64 = 2.0;
    let xnu: f64 = 0.5;
    let _tsw: f64 = 1.0;
    let amax = amaxg * 32.2;
    let mut _beta = xntic;

    let mut rng: rand::rngs::StdRng = match seed {
        Some(s) => rand::rngs::StdRng::seed_from_u64(s),
        None => rand::rngs::StdRng::from_entropy(),
    };
    let normal = Normal::new(0.0, 1.0).unwrap();
    let uniform = Uniform::new(0.0_f64, 1.0_f64);

    let mut array_tf = Vec::new();
    let mut array_rms = Vec::new();

    let mut tf_val = 0.1;
    while tf_val <= 6.0 + 1e-6 {
        let mut z: Vec<f64> = vec![0.0; num_runs];
        let mut z1: f64 = 0.0;

        for jj in 0..num_runs {
            let sum: f64 = uniform.sample(&mut rng);
            let tstart = tf_val * sum;
            let _sum: f64 = uniform.sample(&mut rng);
            let coef = 1.0;
            let _phase = 6.28 * uniform.sample(&mut rng);

            let mut y: f64 = 0.0;
            let mut yd: f64 = 0.0;
            let ts2 = ts * ts;
            let ts3 = ts2 * ts;
            let _ts4 = ts3 * ts;
            let _ts5 = _ts4 * ts;

            // Initialize matrices
            let mut phi = [[0.0; 3]; 3];
            let mut p = [[0.0; 3]; 3];
            let mut qc = [[0.0; 3]; 3];
            let mut idnp = [[0.0; 3]; 3];
            let mut f = [[0.0; 3]; 3];

            let rtm = vc * tf_val;
            let sigpos = rtm * signoise;
            let sign2 = sigpos * sigpos;

            idnp[0][0] = 1.0;
            idnp[1][1] = 1.0;
            idnp[2][2] = 1.0;

            let hmat = [1.0, 0.0, 0.0];

            p[0][0] = sign2;
            p[1][1] = (vm * hedeg / 57.3).powi(2);
            p[2][2] = xntic * xntic;

            let phin: f64;
            if ipoisson == 1 {
                f[0][1] = 1.0;
                f[1][2] = 1.0;
                f[2][2] = -2.0 * xnu;
                phin = 4.0 * xnu * xntic * xntic;
            } else {
                f[0][1] = 1.0;
                f[1][2] = 1.0;
                phin = xntic * xntic / 10.0;
            }

            qc[2][2] = phin;

            // PHI = IDNP + F*TS + 0.5*F*F*TS*TS
            let mut ff = [[0.0; 3]; 3];
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

            // Q = PHI * QC * PHI' * TS
            let mut phiqc = [[0.0; 3]; 3];
            for i in 0..order {
                for j in 0..order {
                    for k in 0..order {
                        phiqc[i][j] += phi[i][k] * qc[k][j];
                    }
                }
            }
            let mut q = [[0.0; 3]; 3];
            for i in 0..order {
                for j in 0..order {
                    for k in 0..order {
                        q[i][j] += phiqc[i][k] * phi[j][k]; // phi' means transposed
                    }
                    q[i][j] *= ts;
                }
            }

            let mut t: f64 = 0.0;
            let h: f64 = 0.001;
            let mut s: f64 = 0.0;
            let mut yh: f64 = 0.0;
            let mut ydh: f64 = 0.0;
            let mut xnth: f64 = 0.0;
            let mut xnc: f64 = 0.0;
            let mut xnl: f64 = 0.0;
            _beta = xntic;
            let mut qfirst = true;
            let sig = 1.0 / (2.0 * xnu).sqrt();
            let xnoise: f64 = normal.sample(&mut rng);
            let mut xntp = if xnoise > 0.0 { _beta } else { -_beta };
            let mut delt = 9999.0;
            let mut tnow = 0.0;

            while t <= tf_val - 0.0001 {
                if qfirst {
                    let xnoise1 = sig * normal.sample(&mut rng);
                    let xnoise2 = sig * normal.sample(&mut rng);
                    delt = xnoise1 * xnoise1 + xnoise2 * xnoise2;
                    qfirst = false;
                    tnow = t;
                }
                if t >= delt + tnow {
                    xntp = -xntp;
                    qfirst = true;
                }

                let yold = y;
                let ydold = yd;
                let xnlold = xnl;

                // First derivative evaluation
                let xntc = if iconstant == 1 {
                    xntp
                } else if iconstant == 2 {
                    if t < tstart {
                        0.0
                    } else {
                        coef * xntic
                    }
                } else if t < tstart {
                    0.0
                } else {
                    xntic * (2.0 * t + _phase).sin()
                };

                let tgo = tf_val - t + 0.00001;
                let _rtm = vc * tgo;
                let xlam = y / (vc * tgo);
                let xnld = (xnc - xnl) / tau;
                let ydd = xntc - xnl;

                // Euler step
                y += h * yd;
                yd += h * ydd;
                xnl += h * xnld;
                t += h;

                // Second derivative for RK2
                let xntc = if iconstant == 1 {
                    xntp
                } else if iconstant == 2 {
                    if t < tstart {
                        0.0
                    } else {
                        coef * xntic
                    }
                } else if t < tstart {
                    0.0
                } else {
                    xntic * (2.0 * t + _phase).sin()
                };

                let tgo = tf_val - t + 0.00001;
                let _rtm = vc * tgo;
                let _xlam = y / (vc * tgo);
                let xnld = (xnc - xnl) / tau;
                let ydd = xntc - xnl;

                // RK2 averaging
                y = 0.5 * (yold + y + h * yd);
                yd = 0.5 * (ydold + yd + h * ydd);
                xnl = 0.5 * (xnlold + xnl + h * xnld);

                s += h;
                if s >= ts - 0.00001 {
                    s = 0.0;
                    let tgo = tf_val - t + 0.000001;
                    let rtm = vc * tgo;
                    let sigpos = rtm * signoise;
                    let sign2 = sigpos * sigpos;
                    let rmat = sign2;

                    // M = PHI * P * PHI' + Q
                    let mut phip = [[0.0; 3]; 3];
                    for i in 0..order {
                        for j in 0..order {
                            for k in 0..order {
                                phip[i][j] += phi[i][k] * p[k][j];
                            }
                        }
                    }
                    let mut m = [[0.0; 3]; 3];
                    for i in 0..order {
                        for j in 0..order {
                            for k in 0..order {
                                m[i][j] += phip[i][k] * phi[j][k];
                            }
                            m[i][j] += q[i][j];
                        }
                    }

                    // K = M * HMAT' / (HMAT * M * HMAT' + RMAT)
                    // Since HMAT = [1 0 0], HMAT*M*HMAT' = M[0][0]
                    let hmh = m[0][0] + rmat;
                    let mut k = [0.0; 3];
                    for i in 0..order {
                        k[i] = m[i][0] / hmh;
                    }

                    // P = (IDNP - K*HMAT) * M
                    let mut kh = [[0.0; 3]; 3];
                    for i in 0..order {
                        kh[i][0] = k[i] * hmat[0];
                    }
                    let mut idkh = [[0.0; 3]; 3];
                    for i in 0..order {
                        for j in 0..order {
                            idkh[i][j] = idnp[i][j] - kh[i][j];
                        }
                    }
                    let mut new_p = [[0.0; 3]; 3];
                    for i in 0..order {
                        for j in 0..order {
                            for kk in 0..order {
                                new_p[i][j] += idkh[i][kk] * m[kk][j];
                            }
                        }
                    }
                    p = new_p;

                    let xlamnoise = signoise * normal.sample(&mut rng);
                    let ystar = rtm * (xlam + xlamnoise);

                    let (yb, ydb, xntb) = project34(t, ts, yh, ydh, xnl, xnth, xnu, ipoisson);

                    let res = ystar - yb;
                    yh = yb + k[0] * res;
                    ydh = ydb + k[1] * res;
                    xnth = xntb + k[2] * res;

                    let x = tgo / tau;
                    let _tautgt = 1.0 / (2.0 * xnu);
                    let zem1h = yh + ydh * tgo - xnl * tau * tau * ((-x).exp() + x - 1.0)
                        + 0.5 * xnth * tgo * tgo;
                    let top = 6.0 * x * x * ((-x).exp() - 1.0 + x);
                    let bot1 = 2.0 * x.powi(3) + 3.0 + 6.0 * x - 6.0 * x * x;
                    let bot2 = -12.0 * x * (-x).exp() - 3.0 * (-2.0 * x).exp();
                    let xnpp = top / (0.0001 + bot1 + bot2);
                    xnc = xnpp * zem1h / tgo.powi(2);
                    if xnc > amax {
                        xnc = amax;
                    }
                    if xnc < -amax {
                        xnc = -amax;
                    }
                }
            }

            z[jj] = y;
            z1 += z[jj];
            let _xmean = z1 / (jj + 1) as f64;
        }

        // Calculate RMS
        let mut z2: f64 = 0.0;
        for jj in 0..num_runs {
            z2 += z[jj] * z[jj];
        }
        let rms = if num_runs > 1 {
            (z2 / (num_runs - 1) as f64).sqrt()
        } else {
            0.0
        };

        array_tf.push(tf_val);
        array_rms.push(rms);

        tf_val += 0.1;
    }

    Results {
        tf: array_tf,
        rms: array_rms,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c34l3_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.tf.clone(),
        results.rms.clone(),
    ])?;

    let plot_file = format!("{}/c34l3_rms.png", output_dir);
    let config = PlotConfig::new("RMS Miss for Various Flight Times")
        .with_labels("Flight Time (S)", "RMS MISS (Ft)");

    let series = vec![
        Series::new(results.tf.clone(), results.rms.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C34L3: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c34l3_runs() {
        let results = run_with_seed(Some(12345));
        assert!(!results.tf.is_empty());
    }
}
