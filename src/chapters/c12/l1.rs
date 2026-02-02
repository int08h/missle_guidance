//! Chapter 12, Lesson 1: Extended Kalman Filter
//!
//! Estimates altitude, velocity, and ballistic coefficient of a
//! reentry vehicle using an Extended Kalman Filter.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;
use rand::prelude::*;
use rand::SeedableRng;
use rand_distr::StandardNormal;

pub struct Results {
    pub time: Vec<f64>,
    pub x: Vec<f64>,
    pub xh: Vec<f64>,
    pub xd: Vec<f64>,
    pub xdh: Vec<f64>,
    pub beta: Vec<f64>,
    pub betah: Vec<f64>,
}

/// Project state forward using simple Euler integration
fn project_c12l1(ts: f64, xp: f64, xdp: f64, beta: f64, hp: f64) -> (f64, f64) {
    let mut t: f64 = 0.0;
    let mut x = xp;
    let mut xd = xdp;
    let h = hp;

    while t <= ts - 0.0001 {
        let xdd = 0.0034 * 32.2 * xd * xd * (-x / 22000.0).exp() / (2.0 * beta) - 32.2;
        xd += h * xdd;
        x += h * xd;
        t += h;
    }

    (x, xd)
}

/// 3x3 matrix operations
fn mat3_mul(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}

fn mat3_transpose(a: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = a[j][i];
        }
    }
    c
}

fn mat3_add(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = a[i][j] + b[i][j];
        }
    }
    c
}

fn mat3_sub(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = a[i][j] - b[i][j];
        }
    }
    c
}

/// Run the C12L1 simulation
pub fn run() -> Results {
    run_with_seed(None)
}

/// Run with optional seed for reproducibility
pub fn run_with_seed(seed: Option<u64>) -> Results {
    let iterm: i32 = 1;
    let signoise: f64 = 25.0;
    let mut x: f64 = 200000.0;
    let mut xd: f64 = -6000.0;
    let beta: f64 = 500.0;
    let mut xh: f64 = 200025.0;
    let mut xdh: f64 = -6150.0;
    let mut betah: f64 = 800.0;
    let ts: f64 = 0.1;
    let tf: f64 = 30.0;
    let phis: f64 = 0.0;
    let h: f64 = 0.001;
    let hp: f64 = 0.001;

    let mut rng: Box<dyn RngCore> = match seed {
        Some(s) => Box::new(rand::rngs::StdRng::seed_from_u64(s)),
        None => Box::new(rand::rngs::StdRng::from_entropy()),
    };

    // Initial covariance
    let mut p = [[0.0f64; 3]; 3];
    p[0][0] = signoise * signoise;
    p[1][1] = 20000.0;
    p[2][2] = 300.0 * 300.0;

    let idnp = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let _hmat = [1.0, 0.0, 0.0];
    let rmat = signoise * signoise;

    let mut t: f64 = 0.0;
    let mut s: f64 = 0.0;

    let mut array_t = Vec::new();
    let mut array_x = Vec::new();
    let mut array_xh = Vec::new();
    let mut array_xd = Vec::new();
    let mut array_xdh = Vec::new();
    let mut array_beta = Vec::new();
    let mut array_betah = Vec::new();

    while t <= tf {
        let xold = x;
        let xdold = xd;

        // Truth model integration
        let xdd = 0.0034 * 32.2 * xd * xd * (-x / 22000.0).exp() / (2.0 * beta) - 32.2;
        x += h * xd;
        xd += h * xdd;
        t += h;

        let xdd = 0.0034 * 32.2 * xd * xd * (-x / 22000.0).exp() / (2.0 * beta) - 32.2;
        x = 0.5 * (xold + x + h * xd);
        xd = 0.5 * (xdold + xd + h * xdd);

        s += h;

        if s >= ts - 0.00001 {
            s = 0.0;

            // Compute F matrix (linearized dynamics)
            let rhoh = 0.0034 * (-xh / 22000.0).exp();
            let f21 = -32.2 * rhoh * xdh * xdh / (44000.0 * betah);
            let f22 = rhoh * 32.2 * xdh / betah;
            let f23 = -rhoh * 32.2 * xdh * xdh / (2.0 * betah * betah);

            let mut f = [[0.0f64; 3]; 3];
            f[0][1] = 1.0;
            f[1][0] = f21;
            f[1][1] = f22;
            f[1][2] = f23;

            // State transition matrix
            let phi = if iterm == 1 {
                let mut phi = idnp;
                for i in 0..3 {
                    for j in 0..3 {
                        phi[i][j] += f[i][j] * ts;
                    }
                }
                phi
            } else {
                // Second order approximation
                let f2 = mat3_mul(&f, &f);
                let mut phi = idnp;
                for i in 0..3 {
                    for j in 0..3 {
                        phi[i][j] += f[i][j] * ts + f2[i][j] * ts * ts / 2.0;
                    }
                }
                phi
            };

            // Process noise
            let mut q = [[0.0f64; 3]; 3];
            q[1][1] = f23 * f23 * phis * ts * ts * ts / 3.0;
            q[1][2] = f23 * phis * ts * ts / 2.0;
            q[2][1] = f23 * phis * ts * ts / 2.0;
            q[2][2] = phis * ts;

            // Predict covariance: M = PHI * P * PHI' + Q
            let phi_t = mat3_transpose(&phi);
            let temp = mat3_mul(&phi, &p);
            let m = mat3_add(&mat3_mul(&temp, &phi_t), &q);

            // Kalman gain: K = M * H' / (H * M * H' + R)
            // Since H = [1, 0, 0], this simplifies
            let hmht_r = m[0][0] + rmat;
            let k = [m[0][0] / hmht_r, m[1][0] / hmht_r, m[2][0] / hmht_r];

            // Update covariance: P = (I - K*H) * M
            let mut kh = [[0.0f64; 3]; 3];
            kh[0][0] = k[0];
            kh[1][0] = k[1];
            kh[2][0] = k[2];
            p = mat3_mul(&mat3_sub(&idnp, &kh), &m);

            // Measurement update
            let xnoise = signoise * rng.sample::<f64, _>(StandardNormal);
            let (xb, xdb) = project_c12l1(ts, xh, xdh, betah, hp);
            let res = x + xnoise - xb;

            xh = xb + k[0] * res;
            xdh = xdb + k[1] * res;
            betah += k[2] * res;

            array_t.push(t);
            array_x.push(x);
            array_xh.push(xh);
            array_xd.push(xd);
            array_xdh.push(xdh);
            array_beta.push(beta);
            array_betah.push(betah);
        }
    }

    Results {
        time: array_t,
        x: array_x,
        xh: array_xh,
        xd: array_xd,
        xdh: array_xdh,
        beta: array_beta,
        betah: array_betah,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c12l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.x.clone(),
        results.xh.clone(),
        results.xd.clone(),
        results.xdh.clone(),
        results.beta.clone(),
        results.betah.clone(),
    ])?;

    let plot_file = format!("{}/c12l1_altitude.png", output_dir);
    let config = PlotConfig::new("EKF - Altitude Estimation")
        .with_labels("Time (Sec)", "Altitude (Ft)");

    let series = vec![
        Series::new(results.time.clone(), results.x.clone())
            .with_color(plotters::prelude::BLUE)
            .with_label("True"),
        Series::new(results.time.clone(), results.xh.clone())
            .with_color(plotters::prelude::RED)
            .with_label("Estimated"),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C12L1: Simulation finished");
    println!("  Data saved to: {}", data_file);
    println!("  Plot saved to: {}", plot_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c12l1_runs() {
        let results = run_with_seed(Some(12345));
        assert!(!results.time.is_empty());
    }

    #[test]
    fn test_c12l1_filter_converges() {
        let results = run_with_seed(Some(12345));
        // Beta estimate should converge towards true value
        let final_betah = *results.betah.last().unwrap();
        assert!((final_betah - 500.0).abs() < 200.0);  // Within 200 of true value
    }
}
