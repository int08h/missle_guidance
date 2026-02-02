//! Chapter 35, Lesson 6: Miss Distance for Various Flight Times
//!
//! Computes miss distance using optimal guidance gains generated
//! from the GENERATEGAINS function.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub tf: Vec<f64>,
    pub y: Vec<f64>,
}

/// Generate optimal guidance gains using Riccati equation
#[allow(clippy::type_complexity)]
fn generate_gains(
    tau: f64,
    w: f64,
    z: f64,
    wz: f64,
    gam: f64,
    tf: f64,
    ts: f64,
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    // Initialize matrices (6x6)
    let mut f = [[0.0; 6]; 6];
    let mut s = [[0.0; 6]; 6];
    let mut g = [0.0; 6];

    g[5] = w * w / tau;

    f[0][1] = 1.0;
    f[1][2] = 1.0;
    f[1][3] = -1.0;
    f[1][5] = 1.0 / (wz * wz);
    f[3][4] = 1.0;
    f[4][5] = 1.0;
    f[5][3] = -w * w / tau;
    f[5][4] = -w * w * (2.0 * z / w + tau) / tau;
    f[5][5] = -w * w * (1.0 / (w * w) + 2.0 * z * tau / w) / tau;

    s[0][0] = 1.0;

    let mut t: f64 = 0.0;
    let h: f64 = 0.0001;
    let mut s1: f64 = 0.0;

    let mut c1_vec = Vec::new();
    let mut c2_vec = Vec::new();
    let mut c3_vec = Vec::new();
    let mut c4_vec = Vec::new();
    let mut c5_vec = Vec::new();
    let mut c6_vec = Vec::new();

    let gaminv = 1.0 / gam;

    while t < tf - 0.0001 {
        s1 += h;
        let sold = s;

        // Compute SF = S * F
        let mut sf = [[0.0; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                for k in 0..6 {
                    sf[i][j] += s[i][k] * f[k][j];
                }
            }
        }

        // Compute GTS = G' * S
        let mut gts = [0.0; 6];
        for j in 0..6 {
            for k in 0..6 {
                gts[j] += g[k] * s[k][j];
            }
        }

        // C = (1/gam) * GTS
        let mut c = [0.0; 6];
        for i in 0..6 {
            c[i] = gaminv * gts[i];
        }

        // SFT = SF'
        let mut sft = [[0.0; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                sft[i][j] = sf[j][i];
            }
        }

        // CTBC = gam * C' * C
        let mut ctbc = [[0.0; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                ctbc[i][j] = gam * c[i] * c[j];
            }
        }

        // SD = SF + SFT - CTBC
        let mut sd = [[0.0; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                sd[i][j] = sf[i][j] + sft[i][j] - ctbc[i][j];
            }
        }

        // Euler step
        for i in 0..6 {
            for j in 0..6 {
                s[i][j] += h * sd[i][j];
            }
        }
        t += h;

        // Recompute for RK2
        let mut sf = [[0.0; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                for k in 0..6 {
                    sf[i][j] += s[i][k] * f[k][j];
                }
            }
        }

        let mut gts = [0.0; 6];
        for j in 0..6 {
            for k in 0..6 {
                gts[j] += g[k] * s[k][j];
            }
        }

        let mut c = [0.0; 6];
        for i in 0..6 {
            c[i] = gaminv * gts[i];
        }

        let mut sft = [[0.0; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                sft[i][j] = sf[j][i];
            }
        }

        let mut ctbc = [[0.0; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                ctbc[i][j] = gam * c[i] * c[j];
            }
        }

        let mut sd = [[0.0; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                sd[i][j] = sf[i][j] + sft[i][j] - ctbc[i][j];
            }
        }

        // RK2 averaging
        for i in 0..6 {
            for j in 0..6 {
                s[i][j] = 0.5 * (sold[i][j] + s[i][j]) + 0.5 * h * sd[i][j];
            }
        }

        if s1 >= ts - 0.0001 {
            s1 = 0.0;
            c1_vec.push(-c[0]);
            c2_vec.push(-c[1]);
            c3_vec.push(-c[2]);
            c4_vec.push(-c[3]);
            c5_vec.push(-c[4]);
            c6_vec.push(-c[5]);
        }
    }

    (c1_vec, c2_vec, c3_vec, c4_vec, c5_vec, c6_vec)
}

/// Run the C35L6 simulation
pub fn run() -> Results {
    let tau: f64 = 1.0;
    let w: f64 = 20.0;
    let z: f64 = 0.7;
    let wz: f64 = 5.0;
    let apn: i32 = 3;
    let xnt: f64 = 16.1;
    let ts: f64 = 0.01;
    let gam: f64 = 0.00001;
    let xlim: f64 = 9999999.0;
    let tfmax: f64 = 10.0;
    let xnp: f64 = 3.0;
    let vc: f64 = 4000.0;
    let _vm: f64 = 3000.0;

    // Generate gains if APN == 3
    let (c1, c2, c3, c4, c5, c6) = if apn == 3 {
        generate_gains(tau, w, z, wz, gam, tfmax, ts)
    } else {
        (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new())
    };

    let mut array_tf = Vec::new();
    let mut array_y = Vec::new();

    let mut tf_val = 0.1;
    while tf_val <= 10.0 + 1e-6 {
        let mut e: f64 = 0.0;
        let mut ed: f64 = 0.0;
        let mut edd: f64 = 0.0;
        let mut t: f64 = 0.0;
        let h: f64 = 0.0001;
        let mut s: f64 = 0.0;
        let mut y: f64 = 0.0;
        let mut yd: f64 = 0.0;
        let mut xnc: f64 = 0.0;
        let _rtm = vc * tf_val;

        while t < tf_val - 0.00001 {
            s += h;
            let eold = e;
            let edold = ed;
            let eddold = edd;
            let yold = y;
            let ydold = yd;

            // First derivative evaluation
            let tgo = tf_val - t + 0.0001;
            let _rtm = vc * tgo;
            let _xlam = y / (vc * tgo);
            let _xncg = xnc / 32.2;
            let eddd = w * w * (xnc - e - (2.0 * z / w + tau) * ed
                - (2.0 * z * tau / w + 1.0 / (w * w)) * edd) / tau;
            let xnl = e - edd / (wz * wz);
            let ydd = xnt - xnl;

            // Euler step
            e += h * ed;
            ed += h * edd;
            edd += h * eddd;
            y += h * yd;
            yd += h * ydd;
            t += h;

            // Second derivative for RK2
            let tgo = tf_val - t + 0.0001;
            let _rtm = vc * tgo;
            let _xlam = y / (vc * tgo);
            let _xncg = xnc / 32.2;
            let eddd = w * w * (xnc - e - (2.0 * z / w + tau) * ed
                - (2.0 * z * tau / w + 1.0 / (w * w)) * edd) / tau;
            let xnl = e - edd / (wz * wz);
            let ydd = xnt - xnl;

            // RK2 averaging
            e = 0.5 * (eold + e + h * ed);
            ed = 0.5 * (edold + ed + h * edd);
            edd = 0.5 * (eddold + edd + h * eddd);
            y = 0.5 * (yold + y + h * yd);
            yd = 0.5 * (ydold + yd + h * ydd);

            if s >= ts - 0.0001 {
                s = 0.0;
                let tgo = tf_val - t + 0.0001;

                xnc = if apn == 0 {
                    xnp * (y + yd * tgo) / (tgo * tgo)
                } else if apn == 1 {
                    xnp * (y + yd * tgo + 0.5 * xnt * tgo * tgo) / (tgo * tgo)
                } else if apn == 2 {
                    let xs = tgo / tau;
                    let top = 6.0 * xs * xs * ((-xs).exp() - 1.0 + xs);
                    let bot1 = 2.0 * xs.powi(3) + 3.0 + 6.0 * xs - 6.0 * xs * xs;
                    let bot2 = -12.0 * xs * (-xs).exp() - 3.0 * (-2.0 * xs).exp();
                    let xnpp = top / (0.0001 + bot1 + bot2);
                    let c1p = xnpp / (tgo * tgo);
                    let c2p = xnpp / tgo;
                    let c3p = 0.5 * xnpp;
                    let c4p = -xnpp * ((-xs).exp() + xs - 1.0) / (xs * xs);
                    c1p * y + c2p * yd + c3p * xnt + c4p * xnl
                } else {
                    // APN == 3: use precomputed gains
                    let jj = (tgo / ts) as usize;
                    let jj = jj.min(c1.len().saturating_sub(1));
                    c1[jj] * y + c2[jj] * yd + c3[jj] * xnt + c4[jj] * e
                        + c5[jj] * ed + c6[jj] * edd
                };

                if xnc > xlim {
                    xnc = xlim;
                } else if xnc < -xlim {
                    xnc = -xlim;
                }
                let _xncg = xnc / 32.2;
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

    let data_file = format!("{}/c35l6_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.tf.clone(),
        results.y.clone(),
    ])?;

    let plot_file = format!("{}/c35l6_miss.png", output_dir);
    let config = PlotConfig::new("Miss Distance for Various Flight Times")
        .with_labels("Flight Time (s)", "Miss (ft)");

    let series = vec![
        Series::new(results.tf.clone(), results.y.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C35L6: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c35l6_runs() {
        let results = run();
        assert!(!results.tf.is_empty());
    }
}
