//! Chapter 35, Lesson 2: Optimal Guidance Gains with Oscillator
//!
//! Computes optimal guidance gains using Riccati equation with
//! oscillating target model.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub c1: Vec<f64>,
    pub c1th: Vec<f64>,
    pub c2: Vec<f64>,
    pub c2th: Vec<f64>,
    pub c3: Vec<f64>,
    pub c3th: Vec<f64>,
    pub c4: Vec<f64>,
    pub c4th: Vec<f64>,
    pub np: Vec<f64>,
    pub xnpp: Vec<f64>,
}

/// Run the C35L2 simulation
pub fn run() -> Results {
    let w: f64 = 1.0;
    let gam: f64 = 0.00001;
    let tf: f64 = 10.0;

    // Initialize matrices (4x4)
    let mut f = [[0.0; 4]; 4];
    let mut s = [[0.0; 4]; 4];
    let g = [0.0, -1.0, 0.0, 0.0];

    f[0][1] = 1.0;
    f[1][2] = 1.0;
    f[2][3] = 1.0;
    f[3][2] = -w * w;

    s[0][0] = 1.0;

    let mut t: f64 = 0.0;
    let h: f64 = 0.0001;
    let mut s1: f64 = 0.0;

    let mut array_t = Vec::new();
    let mut array_c1 = Vec::new();
    let mut array_c1th = Vec::new();
    let mut array_c2 = Vec::new();
    let mut array_c2th = Vec::new();
    let mut array_c3 = Vec::new();
    let mut array_c3th = Vec::new();
    let mut array_c4 = Vec::new();
    let mut array_c4th = Vec::new();
    let mut array_np = Vec::new();
    let mut array_xnpp = Vec::new();

    while t < tf - 0.0001 {
        s1 += h;
        let sold = s;

        // Compute SF = S * F
        let mut sf = [[0.0; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                for k in 0..4 {
                    sf[i][j] += s[i][k] * f[k][j];
                }
            }
        }

        // Compute GTS = G' * S (1x4 row vector)
        let mut gts = [0.0; 4];
        for j in 0..4 {
            for k in 0..4 {
                gts[j] += g[k] * s[k][j];
            }
        }

        // C = (1/gam) * GTS
        let gaminv = 1.0 / gam;
        let c = [
            gaminv * gts[0],
            gaminv * gts[1],
            gaminv * gts[2],
            gaminv * gts[3],
        ];

        // SFT = SF' (transpose of SF)
        let mut sft = [[0.0; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                sft[i][j] = sf[j][i];
            }
        }

        // CTBC = gam * C' * C (outer product scaled by gam)
        let mut ctbc = [[0.0; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                ctbc[i][j] = gam * c[i] * c[j];
            }
        }

        // SD = SF + SFT - CTBC
        let mut sd = [[0.0; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                sd[i][j] = sf[i][j] + sft[i][j] - ctbc[i][j];
            }
        }

        // Euler step: S = S + H * SD
        for i in 0..4 {
            for j in 0..4 {
                s[i][j] += h * sd[i][j];
            }
        }
        t += h;

        // Recompute SD for RK2
        let mut sf = [[0.0; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                for k in 0..4 {
                    sf[i][j] += s[i][k] * f[k][j];
                }
            }
        }

        let mut gts = [0.0; 4];
        for j in 0..4 {
            for k in 0..4 {
                gts[j] += g[k] * s[k][j];
            }
        }

        let c = [
            gaminv * gts[0],
            gaminv * gts[1],
            gaminv * gts[2],
            gaminv * gts[3],
        ];

        let mut sft = [[0.0; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                sft[i][j] = sf[j][i];
            }
        }

        let mut ctbc = [[0.0; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                ctbc[i][j] = gam * c[i] * c[j];
            }
        }

        let mut sd = [[0.0; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                sd[i][j] = sf[i][j] + sft[i][j] - ctbc[i][j];
            }
        }

        // RK2 averaging: S = 0.5*(SOLD + S) + 0.5*H*SD
        for i in 0..4 {
            for j in 0..4 {
                s[i][j] = 0.5 * (sold[i][j] + s[i][j]) + 0.5 * h * sd[i][j];
            }
        }

        if s1 >= 0.009999 {
            s1 = 0.0;
            let c1 = -c[0];
            let c2 = -c[1];
            let c3 = -c[2];
            let c4 = -c[3];
            let np = c2 * t;
            let xnpp = 3.0;

            // Theoretical values
            let c1th = xnpp / (t * t);
            let c2th = xnpp / t;
            let c3th = xnpp * ((1.0 - (w * t).cos()) / w.powi(2)) / t.powi(2);
            let c4th = xnpp * ((w * t - (w * t).sin()) / w.powi(3)) / t.powi(2);

            array_t.push(t);
            array_c1.push(c1);
            array_c1th.push(c1th);
            array_c2.push(c2);
            array_c2th.push(c2th);
            array_c3.push(c3);
            array_c3th.push(c3th);
            array_c4.push(c4);
            array_c4th.push(c4th);
            array_np.push(np);
            array_xnpp.push(xnpp);
        }
    }

    Results {
        time: array_t,
        c1: array_c1,
        c1th: array_c1th,
        c2: array_c2,
        c2th: array_c2th,
        c3: array_c3,
        c3th: array_c3th,
        c4: array_c4,
        c4th: array_c4th,
        np: array_np,
        xnpp: array_xnpp,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c35l2_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.c1.clone(),
        results.c1th.clone(),
        results.c2.clone(),
        results.c2th.clone(),
        results.c3.clone(),
        results.c3th.clone(),
        results.c4.clone(),
        results.c4th.clone(),
        results.np.clone(),
        results.xnpp.clone(),
    ])?;

    let plot_file = format!("{}/c35l2_c3.png", output_dir);
    let config = PlotConfig::new("Optimal Guidance Gain C3")
        .with_labels("Time (s)", "C3");

    let series = vec![
        Series::new(results.time.clone(), results.c3.clone())
            .with_color(plotters::prelude::BLUE)
            .with_label("C3"),
        Series::new(results.time.clone(), results.c3th.clone())
            .with_color(plotters::prelude::RED)
            .with_label("C3 Theory"),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C35L2: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c35l2_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
