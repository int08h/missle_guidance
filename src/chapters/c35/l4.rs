//! Chapter 35, Lesson 4: Riccati with Autopilot Dynamics
//!
//! Computes optimal guidance gains using Riccati equation with
//! third-order autopilot dynamics.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub c4: Vec<f64>,
    pub c5: Vec<f64>,
    pub c6: Vec<f64>,
    pub np: Vec<f64>,
}

/// Run the C35L4 simulation
pub fn run() -> Results {
    let tau: f64 = 1.0;
    let wz: f64 = 5.0;
    let w: f64 = 20.0;
    let z: f64 = 0.7;
    let tf: f64 = 10.0;
    let gam: f64 = 0.00001;

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

    let mut array_t = Vec::new();
    let mut array_c4 = Vec::new();
    let mut array_c5 = Vec::new();
    let mut array_c6 = Vec::new();
    let mut array_np = Vec::new();

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

        // Compute GTS = G' * S (1x6 row vector)
        let mut gts = [0.0; 6];
        for j in 0..6 {
            for k in 0..6 {
                gts[j] += g[k] * s[k][j];
            }
        }

        // C = (1/gam) * GTS
        let gaminv = 1.0 / gam;
        let mut c = [0.0; 6];
        for i in 0..6 {
            c[i] = gaminv * gts[i];
        }

        // SFT = SF' (transpose of SF)
        let mut sft = [[0.0; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                sft[i][j] = sf[j][i];
            }
        }

        // CTBC = gam * C' * C (outer product scaled by gam)
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

        // Euler step: S = S + H * SD
        for i in 0..6 {
            for j in 0..6 {
                s[i][j] += h * sd[i][j];
            }
        }
        t += h;

        // Recompute SD for RK2
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

        // RK2 averaging: S = 0.5*(SOLD + S) + 0.5*H*SD
        for i in 0..6 {
            for j in 0..6 {
                s[i][j] = 0.5 * (sold[i][j] + s[i][j]) + 0.5 * h * sd[i][j];
            }
        }

        if s1 >= 0.009999 {
            s1 = 0.0;
            let _c1 = -c[0];
            let c2 = -c[1];
            let _c3 = -c[2];
            let c4 = -c[3];
            let c5 = -c[4];
            let c6 = -c[5];
            let np = c2 * t;

            array_t.push(t);
            array_c4.push(c4);
            array_c5.push(c5);
            array_c6.push(c6);
            array_np.push(np);
        }
    }

    Results {
        time: array_t,
        c4: array_c4,
        c5: array_c5,
        c6: array_c6,
        np: array_np,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c35l4_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.c4.clone(),
        results.c5.clone(),
        results.c6.clone(),
        results.np.clone(),
    ])?;

    let plot_file = format!("{}/c35l4_np.png", output_dir);
    let config = PlotConfig::new("Navigation Ratio NP")
        .with_labels("Time (s)", "NP");

    let series = vec![
        Series::new(results.time.clone(), results.np.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C35L4: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c35l4_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
