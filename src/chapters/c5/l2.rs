//! Chapter 5, Lesson 2: Covariance Analysis with Matrix Riccati Equation
//!
//! Propagates covariance matrix using RK4 for miss distance analysis.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub sig_y: Vec<f64>,
    pub sig_nl_g: Vec<f64>,
}

/// 4x4 matrix operations
fn mat_add(a: &[[f64; 4]; 4], b: &[[f64; 4]; 4]) -> [[f64; 4]; 4] {
    let mut c = [[0.0; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            c[i][j] = a[i][j] + b[i][j];
        }
    }
    c
}

fn mat_scale(a: &[[f64; 4]; 4], s: f64) -> [[f64; 4]; 4] {
    let mut c = [[0.0; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            c[i][j] = a[i][j] * s;
        }
    }
    c
}

fn mat_mul(a: &[[f64; 4]; 4], b: &[[f64; 4]; 4]) -> [[f64; 4]; 4] {
    let mut c = [[0.0; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            for k in 0..4 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}

fn mat_transpose(a: &[[f64; 4]; 4]) -> [[f64; 4]; 4] {
    let mut c = [[0.0; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            c[i][j] = a[j][i];
        }
    }
    c
}

/// Run the C5L2 simulation - Covariance analysis
pub fn run() -> Results {
    let xnp: f64 = 3.0;       // Navigation ratio
    let tau: f64 = 1.0;       // Time constant
    let xnt: f64 = 96.6;      // Target acceleration
    let vc: f64 = 4000.0;     // Closing velocity
    let tf: f64 = 10.0;       // Flight time
    let h: f64 = 0.01;        // Time step

    let mut t: f64 = 0.0;
    let mut s: f64 = 0.0;
    let phis = xnt * xnt / tf;

    // Initialize matrices
    let mut f = [[0.0f64; 4]; 4];
    let mut q = [[0.0f64; 4]; 4];
    let mut x = [[0.0f64; 4]; 4];

    f[0][1] = 1.0;
    f[1][2] = 1.0;
    q[2][2] = phis;

    let mut array_t = Vec::new();
    let mut array_sig_y = Vec::new();
    let mut array_sig_nl_g = Vec::new();

    while t < tf - 0.0001 {
        s += h;
        let xold = x;

        // RK4 for matrix differential equation
        // Must update TGO at each substep like MATLAB does
        let compute_xd = |x: &[[f64; 4]; 4], tgo: f64| -> [[f64; 4]; 4] {
            let mut f_local = f;
            f_local[1][0] = -xnp / (tau * tgo);
            f_local[1][3] = xnp * vc / tau;
            f_local[3][0] = 1.0 / (tau * vc * tgo);
            f_local[3][3] = -1.0 / tau;

            let fx = mat_mul(&f_local, x);
            let fxt = mat_transpose(&fx);
            mat_add(&mat_add(&fx, &fxt), &q)
        };

        // k0: evaluated at t
        let tgo0 = tf - t + 0.00001;
        let k0 = compute_xd(&x, tgo0);

        // k1: evaluated at t + 0.5*h
        let tgo1 = tf - (t + 0.5 * h) + 0.00001;
        let x_temp = mat_add(&xold, &mat_scale(&k0, 0.5 * h));
        let k1 = compute_xd(&x_temp, tgo1);

        // k2: evaluated at t + 0.5*h (same TGO as k1)
        let x_temp = mat_add(&xold, &mat_scale(&k1, 0.5 * h));
        let k2 = compute_xd(&x_temp, tgo1);

        // k3: evaluated at t + h
        let tgo3 = tf - (t + h) + 0.00001;
        let x_temp = mat_add(&xold, &mat_scale(&k2, h));
        let k3 = compute_xd(&x_temp, tgo3);

        // Update
        t += h;
        let k_sum = mat_add(
            &mat_add(&k0, &mat_scale(&mat_add(&k1, &k2), 2.0)),
            &k3
        );
        x = mat_add(&xold, &mat_scale(&k_sum, h / 6.0));

        // Sample output
        if s >= 0.09999 {
            s = 0.0;
            let tgo = tf - t + 0.00001;

            // Compute A matrix for acceleration variance
            let mut a = [0.0f64; 4];
            a[0] = xnp / (tau * tgo);
            a[3] = -xnp * vc / tau;

            // AXAT = A * X * A'
            let mut ax = [0.0f64; 4];
            for j in 0..4 {
                for i in 0..4 {
                    ax[j] += a[i] * x[i][j];
                }
            }
            let mut axat = 0.0f64;
            for i in 0..4 {
                axat += ax[i] * a[i];
            }

            let sig_y = x[0][0].sqrt();
            let sig_nl = axat.sqrt() / 32.2;

            array_t.push(t);
            array_sig_y.push(sig_y);
            array_sig_nl_g.push(sig_nl);
        }
    }

    println!("Final SIGY: {:.4}", array_sig_y.last().unwrap_or(&0.0));

    Results {
        time: array_t,
        sig_y: array_sig_y,
        sig_nl_g: array_sig_nl_g,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c5l2_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.sig_y.clone(),
        results.sig_nl_g.clone(),
    ])?;

    let plot_file1 = format!("{}/c5l2_sigy.png", output_dir);
    let config = PlotConfig::new("Covariance Analysis")
        .with_labels("Time (s)", "Std Dev of Relative Position (Ft)");
    let series = vec![Series::new(results.time.clone(), results.sig_y.clone())
        .with_color(plotters::prelude::BLUE)];
    line_plot(&plot_file1, &config, &series).ok();

    let plot_file2 = format!("{}/c5l2_signl.png", output_dir);
    let config = PlotConfig::new("Covariance Analysis")
        .with_labels("Time (s)", "Std Dev of Acceleration (G)")
        .with_y_range(0.0, 20.0);
    let series = vec![Series::new(results.time.clone(), results.sig_nl_g.clone())
        .with_color(plotters::prelude::BLUE)];
    line_plot(&plot_file2, &config, &series).ok();

    println!("C5L2: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c5l2_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
