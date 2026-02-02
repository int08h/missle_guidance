//! Chapter 29, Lesson 5: Adjoint Target Weave FFT Analysis
//!
//! Uses FFT to compute the frequency response (PZ) of the adjoint
//! model for target weaving.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

// Use MATLAB's PI value for exact matching
const PI: f64 = 3.1416;

pub struct Results {
    pub w: Vec<f64>,
    pub pz: Vec<f64>,
}

/// Run the C29L5 simulation
pub fn run() -> Results {
    let xnt: f64 = 96.6;
    let xnp: f64 = 4.0;
    let tau: f64 = 1.0;
    let vm: f64 = 3000.0;
    let hedeg: f64 = -20.0;

    let tf: f64 = 25.6;
    let ts: f64 = 0.1;
    let xnpoint = tf / ts;
    let npoint = xnpoint.round() as usize;
    let fs = 1.0 / ts;

    let t: f64 = 0.0;
    let mut s: f64 = 0.0;
    let mut tp = t + 0.00001;
    let mut x1: f64 = 0.0;
    let mut x2: f64 = 0.0;
    let mut x3: f64 = 1.0;
    let mut x4: f64 = 0.0;
    let h: f64 = 0.01;
    let _he = hedeg / 57.3;

    let mut x_arr = vec![0.0; npoint];
    let mut j: usize = 0;

    while tp <= tf - 0.0001 {
        let x1old = x1;
        let x2old = x2;
        let x3old = x3;
        let x4old = x4;

        // First derivative evaluation
        let x1d = x2;
        let x2d = x3;
        let y1 = (x4 - x2) / tau;
        let tgo = tp + 0.00001;
        let x3d = xnp * y1 / tgo;
        let x4d = -y1;

        // Euler step
        x1 += h * x1d;
        x2 += h * x2d;
        x3 += h * x3d;
        x4 += h * x4d;
        tp += h;

        // Second derivative evaluation
        let x1d = x2;
        let x2d = x3;
        let y1 = (x4 - x2) / tau;
        let tgo = tp + 0.00001;
        let x3d = xnp * y1 / tgo;
        let x4d = -y1;

        // RK2 averaging
        x1 = (x1old + x1) / 2.0 + 0.5 * h * x1d;
        x2 = (x2old + x2) / 2.0 + 0.5 * h * x2d;
        x3 = (x3old + x3) / 2.0 + 0.5 * h * x3d;
        x4 = (x4old + x4) / 2.0 + 0.5 * h * x4d;

        s += h;
        if s >= ts - 0.00001 {
            s = 0.0;
            let _xmnt = xnt * x1;
            let _xmhe = -vm * _he * x2;
            if j < npoint {
                x_arr[j] = x2;
                j += 1;
            }
        }
    }

    // DFT calculation
    let mut xr = vec![0.0; npoint];
    let mut xi = vec![0.0; npoint];

    for k in 0..npoint {
        for (n, &x_val) in x_arr.iter().enumerate() {
            let ag = 2.0 * PI * (k + 1) as f64 * (n + 1) as f64 / npoint as f64;
            xr[k] += x_val * ag.cos() / npoint as f64;
            xi[k] -= x_val * ag.sin() / npoint as f64;
        }
    }

    // Compute magnitude and PZ
    let imax = npoint / 2;
    let mut array_w = Vec::with_capacity(imax + 1);
    let mut array_pz = Vec::with_capacity(imax + 1);

    for i in 0..=imax {
        let f = fs * (i + 1) as f64 / npoint as f64;
        let xmag = (xr[i] * xr[i] + xi[i] * xi[i]).sqrt();
        let w = 2.0 * PI * f;
        let pz = xmag * xnt * xnpoint / fs;
        array_w.push(w);
        array_pz.push(pz);
    }

    Results {
        w: array_w,
        pz: array_pz,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c29l5_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.w.clone(),
        results.pz.clone(),
    ])?;

    let plot_file = format!("{}/c29l5_pz.png", output_dir);
    let config = PlotConfig::new("PZ Frequency Response")
        .with_labels("W (r/s)", "PZ");

    let series = vec![
        Series::new(results.w.clone(), results.pz.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C29L5: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c29l5_runs() {
        let results = run();
        assert!(!results.w.is_empty());
    }
}
