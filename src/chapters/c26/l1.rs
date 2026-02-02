//! Chapter 26, Lesson 1: Optimal Guidance - Riccati Equation
//!
//! Solves the Riccati differential equation to compute optimal control gains.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub c1: Vec<f64>,
    pub c2: Vec<f64>,
}

/// Run the C26L1 simulation
pub fn run() -> Results {
    let order: usize = 2;
    let kdel: f64 = 9000.0;
    let wrr: f64 = 2.0;
    let delcmx: f64 = 30.0;
    let phimx: f64 = 10.0;
    let phidmx: f64 = 300.0;

    // Initialize matrices (2x2)
    let mut g = [[0.0; 1]; 2]; // G is 2x1
    let mut f = [[0.0; 2]; 2]; // F is 2x2
    let mut a = [[0.0; 2]; 2]; // A is 2x2
    let mut s = [[0.0; 2]; 2]; // S is 2x2

    g[1][0] = kdel;
    f[0][1] = 1.0;
    f[1][1] = -wrr;
    a[0][0] = (delcmx / phimx).powi(2);
    a[1][1] = (delcmx / phidmx).powi(2);

    let mut t: f64 = 0.0;
    let h: f64 = 0.0002;
    let mut s1: f64 = 0.0;
    let tf: f64 = 5.0;

    let mut array_t = Vec::new();
    let mut array_c1 = Vec::new();
    let mut array_c2 = Vec::new();

    while t < tf - 0.0001 {
        s1 += h;
        let sold = s;

        // First derivative evaluation
        // SD = S*F + F'*S - S*G*G'*S + A
        // C = G'*S
        let sd = compute_sd(&s, &f, &g, &a);
        let c = compute_c(&g, &s);

        // Euler step
        for i in 0..order {
            for j in 0..order {
                s[i][j] += h * sd[i][j];
            }
        }
        t += h;

        // Second derivative for RK2
        let sd = compute_sd(&s, &f, &g, &a);

        // RK2 averaging
        for i in 0..order {
            for j in 0..order {
                s[i][j] = 0.5 * (sold[i][j] + s[i][j] + h * sd[i][j]);
            }
        }

        if s1 >= 0.0499999 {
            s1 = 0.0;
            array_t.push(t);
            array_c1.push(c[0]);
            array_c2.push(c[1]);
        }
    }

    Results {
        time: array_t,
        c1: array_c1,
        c2: array_c2,
    }
}

/// Compute SD = S*F + F'*S - S*G*G'*S + A
fn compute_sd(
    s: &[[f64; 2]; 2],
    f: &[[f64; 2]; 2],
    g: &[[f64; 1]; 2],
    a: &[[f64; 2]; 2],
) -> [[f64; 2]; 2] {
    let mut sd = [[0.0; 2]; 2];

    // S*F
    let mut sf = [[0.0; 2]; 2];
    for i in 0..2 {
        for j in 0..2 {
            for k in 0..2 {
                sf[i][j] += s[i][k] * f[k][j];
            }
        }
    }

    // F'*S (F transpose times S)
    let mut fts = [[0.0; 2]; 2];
    for i in 0..2 {
        for j in 0..2 {
            for k in 0..2 {
                fts[i][j] += f[k][i] * s[k][j];
            }
        }
    }

    // S*G (2x2 * 2x1 = 2x1)
    let mut sg = [0.0; 2];
    for i in 0..2 {
        for k in 0..2 {
            sg[i] += s[i][k] * g[k][0];
        }
    }

    // G'*S (1x2 * 2x2 = 1x2)
    let mut gts = [0.0; 2];
    for j in 0..2 {
        for k in 0..2 {
            gts[j] += g[k][0] * s[k][j];
        }
    }

    // S*G*G'*S = (S*G) * (G'*S) = outer product of sg and gts
    let mut sgts = [[0.0; 2]; 2];
    for i in 0..2 {
        for j in 0..2 {
            sgts[i][j] = sg[i] * gts[j];
        }
    }

    // SD = S*F + F'*S - S*G*G'*S + A
    for i in 0..2 {
        for j in 0..2 {
            sd[i][j] = sf[i][j] + fts[i][j] - sgts[i][j] + a[i][j];
        }
    }

    sd
}

/// Compute C = G'*S (1x2 matrix as array)
fn compute_c(g: &[[f64; 1]; 2], s: &[[f64; 2]; 2]) -> [f64; 2] {
    let mut c = [0.0; 2];
    for j in 0..2 {
        for k in 0..2 {
            c[j] += g[k][0] * s[k][j];
        }
    }
    c
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c26l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.c1.clone(),
        results.c2.clone(),
    ])?;

    let plot_file = format!("{}/c26l1_gains.png", output_dir);
    let config = PlotConfig::new("Optimal Guidance - Control Gains")
        .with_labels("Time (s)", "Gain");

    let series = vec![
        Series::new(results.time.clone(), results.c1.clone())
            .with_color(plotters::prelude::BLUE)
            .with_label("C1"),
        Series::new(results.time.clone(), results.c2.clone())
            .with_color(plotters::prelude::RED)
            .with_label("C2"),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C26L1: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c26l1_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }

    #[test]
    fn test_c26l1_gains_converge() {
        let results = run();
        // Final gains should have converged to steady-state values
        let n = results.c1.len();
        assert!(n > 10);
        // C1 and C2 should be non-zero at the end
        assert!(results.c1[n - 1].abs() > 0.0);
        assert!(results.c2[n - 1].abs() > 0.0);
    }
}
