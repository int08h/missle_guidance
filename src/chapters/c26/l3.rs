//! Chapter 26, Lesson 3: Optimal Guidance - 4th Order Riccati Equation
//!
//! Solves the 4th order Riccati differential equation to compute optimal control gains
//! with actuator dynamics.
//!
//! MATLAB BUG: Line 41 uses `C=GT*S` but `GT` is never defined. It should be `G'*S`
//! (G transpose times S). This implementation correctly computes G'*S.
//! Output: 5 columns [T, C1, C2, C3, C4]

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub c1: Vec<f64>,
    pub c2: Vec<f64>,
    pub c3: Vec<f64>,
    pub c4: Vec<f64>,
}

const ORDER: usize = 4;

/// Run the C26L3 simulation
pub fn run() -> Results {
    let wact: f64 = 100.0;
    let zact: f64 = 0.65;
    let kdel: f64 = 9000.0;
    let wr: f64 = 2.0;
    let delcmx: f64 = 30.0;
    let phimx: f64 = 10.0;
    let phidmx: f64 = 300.0;

    // Initialize matrices
    let mut g = [[0.0; 1]; ORDER];
    let mut a = [[0.0; ORDER]; ORDER];
    let mut f = [[0.0; ORDER]; ORDER];
    let mut s = [[0.0; ORDER]; ORDER];

    g[3][0] = wact * wact;
    f[0][1] = 1.0;
    f[1][1] = -wr;
    f[1][2] = kdel;
    f[2][3] = 1.0;
    f[3][2] = -wact * wact;
    f[3][3] = -2.0 * zact * wact;
    a[0][0] = (delcmx / phimx).powi(2);
    a[1][1] = (delcmx / phidmx).powi(2);

    let mut t: f64 = 0.0;
    let h: f64 = 0.0002;
    let mut s1: f64 = 0.0;
    let tf: f64 = 5.0;

    let mut array_t = Vec::new();
    let mut array_c1 = Vec::new();
    let mut array_c2 = Vec::new();
    let mut array_c3 = Vec::new();
    let mut array_c4 = Vec::new();

    while t < tf - 0.0001 {
        s1 += h;
        let sold = s;

        // First derivative evaluation
        let sd = compute_sd(&s, &f, &g, &a);
        let c = compute_c(&g, &s);

        // Euler step
        for i in 0..ORDER {
            for j in 0..ORDER {
                s[i][j] += h * sd[i][j];
            }
        }
        t += h;

        // Second derivative for RK2
        let sd = compute_sd(&s, &f, &g, &a);

        // RK2 averaging
        for i in 0..ORDER {
            for j in 0..ORDER {
                s[i][j] = 0.5 * (sold[i][j] + s[i][j] + h * sd[i][j]);
            }
        }

        if s1 >= 0.0499999 {
            s1 = 0.0;
            array_t.push(t);
            array_c1.push(c[0]);
            array_c2.push(c[1]);
            array_c3.push(c[2]);
            array_c4.push(c[3]);
        }
    }

    Results {
        time: array_t,
        c1: array_c1,
        c2: array_c2,
        c3: array_c3,
        c4: array_c4,
    }
}

/// Compute SD = S*F + F'*S - S*G*G'*S + A
fn compute_sd(
    s: &[[f64; ORDER]; ORDER],
    f: &[[f64; ORDER]; ORDER],
    g: &[[f64; 1]; ORDER],
    a: &[[f64; ORDER]; ORDER],
) -> [[f64; ORDER]; ORDER] {
    let mut sd = [[0.0; ORDER]; ORDER];

    // S*F
    let mut sf = [[0.0; ORDER]; ORDER];
    for i in 0..ORDER {
        for j in 0..ORDER {
            for k in 0..ORDER {
                sf[i][j] += s[i][k] * f[k][j];
            }
        }
    }

    // F'*S (F transpose times S)
    let mut fts = [[0.0; ORDER]; ORDER];
    for i in 0..ORDER {
        for j in 0..ORDER {
            for k in 0..ORDER {
                fts[i][j] += f[k][i] * s[k][j];
            }
        }
    }

    // S*G (ORDER x 1)
    let mut sg = [0.0; ORDER];
    for i in 0..ORDER {
        for k in 0..ORDER {
            sg[i] += s[i][k] * g[k][0];
        }
    }

    // G'*S (1 x ORDER)
    let mut gts = [0.0; ORDER];
    for j in 0..ORDER {
        for k in 0..ORDER {
            gts[j] += g[k][0] * s[k][j];
        }
    }

    // S*G*G'*S = outer product of sg and gts
    let mut sgts = [[0.0; ORDER]; ORDER];
    for i in 0..ORDER {
        for j in 0..ORDER {
            sgts[i][j] = sg[i] * gts[j];
        }
    }

    // SD = S*F + F'*S - S*G*G'*S + A
    for i in 0..ORDER {
        for j in 0..ORDER {
            sd[i][j] = sf[i][j] + fts[i][j] - sgts[i][j] + a[i][j];
        }
    }

    sd
}

/// Compute C = G'*S (1 x ORDER matrix as array)
fn compute_c(g: &[[f64; 1]; ORDER], s: &[[f64; ORDER]; ORDER]) -> [f64; ORDER] {
    let mut c = [0.0; ORDER];
    for j in 0..ORDER {
        for k in 0..ORDER {
            c[j] += g[k][0] * s[k][j];
        }
    }
    c
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c26l3_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.c1.clone(),
        results.c2.clone(),
        results.c3.clone(),
        results.c4.clone(),
    ])?;

    let plot_file = format!("{}/c26l3_gains.png", output_dir);
    let config = PlotConfig::new("Optimal Guidance - 4th Order Control Gains")
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

    println!("C26L3: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c26l3_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
