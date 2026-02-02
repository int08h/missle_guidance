//! Chapter 6, Lesson 4: Fifth-Order Binomial Adjoint Noise Analysis
//!
//! Adjoint analysis for noise miss using fifth-order binomial system.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub xmfn: Vec<f64>,
    pub xmrn: Vec<f64>,
    pub xmrna: Vec<f64>,
    pub xmgl: Vec<f64>,
}

/// Run the C6L4 simulation
pub fn run() -> Results {
    let xnp: f64 = 3.0;
    let tau: f64 = 1.0;
    let tf: f64 = 10.0;
    let vc: f64 = 1.0;
    let phifn: f64 = 1.0;
    let phirn: f64 = 1.0;
    let phirna: f64 = 1.0;
    let phigl: f64 = 1.0;
    let ra: f64 = 1.0;
    let h: f64 = 0.01;

    let mut tp: f64 = 0.00001;
    let mut s: f64 = 0.0;

    // State variables (no x1 in this simulation)
    let mut x2: f64 = 0.0;
    let mut x3: f64 = 1.0;
    let mut x4: f64 = 0.0;
    let mut x5: f64 = 0.0;
    let mut x6: f64 = 0.0;
    let mut x7: f64 = 0.0;
    let mut x8: f64 = 0.0;
    let mut x9: f64 = 0.0;
    let mut x10: f64 = 0.0;
    let mut x11: f64 = 0.0;
    let mut x12: f64 = 0.0;

    let mut array_tp = Vec::new();
    let mut array_xmfn = Vec::new();
    let mut array_xmrn = Vec::new();
    let mut array_xmrna = Vec::new();
    let mut array_xmgl = Vec::new();

    while tp <= (tf - 1e-5) {
        let x2_old = x2;
        let x3_old = x3;
        let x4_old = x4;
        let x5_old = x5;
        let x6_old = x6;
        let x7_old = x7;
        let x8_old = x8;
        let x9_old = x9;
        let x10_old = x10;
        let x11_old = x11;
        let x12_old = x12;

        // First derivative evaluation
        let x2d = x3;
        let y1 = 5.0 * (5.0 * x5 / tau + x4) / tau;
        let tgo = tp + 0.00001;
        let x3d = y1 / (vc * tgo);
        let x4d = -y1;
        let x5d = -5.0 * x5 / tau + 5.0 * x6 * xnp * vc / tau;
        let x6d = -5.0 * x6 / tau + 5.0 * x7 / tau;
        let x7d = -5.0 * x7 / tau + 5.0 * x8 / tau;
        let x8d = -5.0 * x8 / tau - x2;
        let x9d = y1 * y1;
        let x10d = (y1 * vc * tgo / ra).powi(2);
        let x11d = (y1 * (vc * tgo / ra).powi(2)).powi(2);
        let x12d = (y1 / (vc * tgo)).powi(2);

        // Euler step
        x2 += h * x2d;
        x3 += h * x3d;
        x4 += h * x4d;
        x5 += h * x5d;
        x6 += h * x6d;
        x7 += h * x7d;
        x8 += h * x8d;
        x9 += h * x9d;
        x10 += h * x10d;
        x11 += h * x11d;
        x12 += h * x12d;
        tp += h;

        // Second derivative evaluation
        let x2d = x3;
        let y1 = 5.0 * (5.0 * x5 / tau + x4) / tau;
        let tgo = tp + 0.00001;
        let x3d = y1 / (vc * tgo);
        let x4d = -y1;
        let x5d = -5.0 * x5 / tau + 5.0 * x6 * xnp * vc / tau;
        let x6d = -5.0 * x6 / tau + 5.0 * x7 / tau;
        let x7d = -5.0 * x7 / tau + 5.0 * x8 / tau;
        let x8d = -5.0 * x8 / tau - x2;
        let x9d = y1 * y1;
        let x10d = (y1 * vc * tgo / ra).powi(2);
        let x11d = (y1 * (vc * tgo / ra).powi(2)).powi(2);
        let x12d = (y1 / (vc * tgo)).powi(2);

        // RK2 averaging
        x2 = 0.5 * (x2_old + x2 + h * x2d);
        x3 = 0.5 * (x3_old + x3 + h * x3d);
        x4 = 0.5 * (x4_old + x4 + h * x4d);
        x5 = 0.5 * (x5_old + x5 + h * x5d);
        x6 = 0.5 * (x6_old + x6 + h * x6d);
        x7 = 0.5 * (x7_old + x7 + h * x7d);
        x8 = 0.5 * (x8_old + x8 + h * x8d);
        x9 = 0.5 * (x9_old + x9 + h * x9d);
        x10 = 0.5 * (x10_old + x10 + h * x10d);
        x11 = 0.5 * (x11_old + x11 + h * x11d);
        x12 = 0.5 * (x12_old + x12 + h * x12d);

        s += h;
        if s >= 0.0999 {
            s = 0.0;
            let xmfn = (x9 * phifn).sqrt();
            let xmrn = (x10 * phirn).sqrt();
            let xmrna = (x11 * phirna).sqrt();
            let xmgl = (x12 * phigl).sqrt();

            array_tp.push(tp);
            array_xmfn.push(xmfn);
            array_xmrn.push(xmrn);
            array_xmrna.push(xmrna);
            array_xmgl.push(xmgl);
        }
    }

    Results {
        time: array_tp,
        xmfn: array_xmfn,
        xmrn: array_xmrn,
        xmrna: array_xmrna,
        xmgl: array_xmgl,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c6l4_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.xmfn.clone(),
        results.xmrn.clone(),
        results.xmrna.clone(),
        results.xmgl.clone(),
    ])?;

    let plot_file = format!("{}/c6l4_fading.png", output_dir);
    let config = PlotConfig::new("Normalized Fading Noise Miss")
        .with_labels("Normalized Flight Time (Sec)", "Normalized Miss");
    let series = vec![
        Series::new(results.time.clone(), results.xmfn.clone())
            .with_color(plotters::prelude::BLUE),
    ];
    line_plot(&plot_file, &config, &series).ok();

    let plot_file2 = format!("{}/c6l4_semiactive.png", output_dir);
    let config2 = PlotConfig::new("Semiactive Noise Miss")
        .with_labels("Normalized Flight Time (Sec)", "Normalized Miss");
    let series2 = vec![
        Series::new(results.time.clone(), results.xmrn.clone())
            .with_color(plotters::prelude::BLUE),
    ];
    line_plot(&plot_file2, &config2, &series2).ok();

    let plot_file3 = format!("{}/c6l4_active.png", output_dir);
    let config3 = PlotConfig::new("Active Noise Miss")
        .with_labels("Normalized Flight Time (Sec)", "Normalized Miss");
    let series3 = vec![
        Series::new(results.time.clone(), results.xmrna.clone())
            .with_color(plotters::prelude::BLUE),
    ];
    line_plot(&plot_file3, &config3, &series3).ok();

    let plot_file4 = format!("{}/c6l4_glint.png", output_dir);
    let config4 = PlotConfig::new("Glint Noise Miss")
        .with_labels("Normalized Flight Time (Sec)", "Normalized Miss");
    let series4 = vec![
        Series::new(results.time.clone(), results.xmgl.clone())
            .with_color(plotters::prelude::BLUE),
    ];
    line_plot(&plot_file4, &config4, &series4).ok();

    println!("C6L4: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c6l4_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
