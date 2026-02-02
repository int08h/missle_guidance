//! Chapter 4, Lesson 7: Adjoint Noise Analysis
//!
//! Determines miss distance standard deviation due to various noise sources
//! (glint, range-independent noise, semiactive/active range-dependent noise)
//! using adjoint equations.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

/// Simulation results
pub struct Results {
    pub tp: Vec<f64>,
    pub xmgl: Vec<f64>,
    pub xmgl_th: Vec<f64>,
    pub xmfn: Vec<f64>,
    pub xmfn_th: Vec<f64>,
    pub xmrn: Vec<f64>,
    pub xmrn_th: Vec<f64>,
    pub xmrna: Vec<f64>,
    pub xmrna_th: Vec<f64>,
}

/// Run the C4L7 simulation
pub fn run() -> Results {
    let xnp: f64 = 3.0;
    let tau: f64 = 0.5;
    let tf: f64 = 10.0;
    let vc: f64 = 4000.0;
    let ts: f64 = 0.1;
    let ra: f64 = 30000.0;

    // Noise parameters
    let siggl: f64 = 10.0;
    let sigfn: f64 = 0.002;
    let sigrn: f64 = 0.02;
    let sigrna: f64 = 0.02;

    let phigl: f64 = siggl * siggl * ts;
    let phifn: f64 = sigfn * sigfn * ts;
    let phirn: f64 = sigrn * sigrn * ts;
    let phirna: f64 = sigrna * sigrna * ts;

    let h: f64 = 0.01;

    let mut tp = 0.00001;
    let mut s = 0.0;

    let mut x1 = 0.0;
    let mut x2 = 0.0;
    let mut x3 = 1.0;
    let mut x4 = 0.0;
    let mut x5 = 0.0;
    let mut x6 = 0.0;
    let mut x7 = 0.0;
    let mut x8 = 0.0;

    let mut array_tp = Vec::new();
    let mut array_xmgl = Vec::new();
    let mut array_xmgl_th = Vec::new();
    let mut array_xmfn = Vec::new();
    let mut array_xmfn_th = Vec::new();
    let mut array_xmrn = Vec::new();
    let mut array_xmrn_th = Vec::new();
    let mut array_xmrna = Vec::new();
    let mut array_xmrna_th = Vec::new();

    while tp <= (tf - 1e-5) {
        s += h;

        let x1_old = x1;
        let x2_old = x2;
        let x3_old = x3;
        let x4_old = x4;
        let x5_old = x5;
        let x6_old = x6;
        let x7_old = x7;
        let x8_old = x8;

        // First step of RK2
        let x1d = x2;
        let x2d = x3;
        let y1 = (x4 - x2) / tau;
        let tgo = tp + 0.00001;
        let x3d = xnp * y1 / tgo;
        let x4d = -y1;
        let x5d = x3d * x3d;
        let x6d = (xnp * vc * y1).powi(2);
        let x7d = (xnp * vc * vc * y1 * tgo / ra).powi(2);
        let x8d = (xnp * vc * vc * vc * y1 * tgo * tgo / (ra * ra)).powi(2);

        x1 += h * x1d;
        x2 += h * x2d;
        x3 += h * x3d;
        x4 += h * x4d;
        x5 += h * x5d;
        x6 += h * x6d;
        x7 += h * x7d;
        x8 += h * x8d;
        tp += h;

        // Second step of RK2
        let x1d = x2;
        let x2d = x3;
        let y1 = (x4 - x2) / tau;
        let tgo = tp + 0.00001;
        let x3d = xnp * y1 / tgo;
        let x4d = -y1;
        let x5d = x3d * x3d;
        let x6d = (xnp * vc * y1).powi(2);
        let x7d = (xnp * vc * vc * y1 * tgo / ra).powi(2);
        let x8d = (xnp * vc * vc * vc * y1 * tgo * tgo / (ra * ra)).powi(2);

        x1 = (x1_old + x1) / 2.0 + 0.5 * h * x1d;
        x2 = (x2_old + x2) / 2.0 + 0.5 * h * x2d;
        x3 = (x3_old + x3) / 2.0 + 0.5 * h * x3d;
        x4 = (x4_old + x4) / 2.0 + 0.5 * h * x4d;
        x5 = (x5_old + x5) / 2.0 + 0.5 * h * x5d;
        x6 = (x6_old + x6) / 2.0 + 0.5 * h * x6d;
        x7 = (x7_old + x7) / 2.0 + 0.5 * h * x7d;
        x8 = (x8_old + x8) / 2.0 + 0.5 * h * x8d;

        s += h;
        if s >= 0.0999 {
            s = 0.0;
            let xmgl = (phigl * x5).sqrt();
            let xmfn = (phifn * x6).sqrt();
            let xmrn = (phirn * x7).sqrt();
            let xmrna = (phirna * x8).sqrt();

            // Theoretical values
            let xmgl_th = 1.44 * (phigl / tau).sqrt();
            let xmfn_th = 0.532 * vc * (tau * phifn).sqrt();
            let xmrn_th = 1.06 * vc * vc * tau.powf(1.5) * phirn.sqrt() / ra;
            let xmrna_th = 4.66 * vc * vc * vc * tau.powf(2.5) * phirna.sqrt() / (ra * ra);

            array_tp.push(tp);
            array_xmgl.push(xmgl);
            array_xmgl_th.push(xmgl_th);
            array_xmfn.push(xmfn);
            array_xmfn_th.push(xmfn_th);
            array_xmrn.push(xmrn);
            array_xmrn_th.push(xmrn_th);
            array_xmrna.push(xmrna);
            array_xmrna_th.push(xmrna_th);
        }
    }

    Results {
        tp: array_tp,
        xmgl: array_xmgl,
        xmgl_th: array_xmgl_th,
        xmfn: array_xmfn,
        xmfn_th: array_xmfn_th,
        xmrn: array_xmrn,
        xmrn_th: array_xmrn_th,
        xmrna: array_xmrna,
        xmrna_th: array_xmrna_th,
    }
}

/// Run and save results to file
pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    // Save data file
    let data_file = format!("{}/c4l7_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.tp.clone(),
        results.xmgl.clone(),
        results.xmgl_th.clone(),
        results.xmfn.clone(),
        results.xmfn_th.clone(),
        results.xmrn.clone(),
        results.xmrn_th.clone(),
        results.xmrna.clone(),
        results.xmrna_th.clone(),
    ])?;

    // Create plots
    let plot_file1 = format!("{}/c4l7_glint.png", output_dir);
    let config1 = PlotConfig::new("Glint Miss Standard Deviation")
        .with_labels("Flight Time (Sec)", "Glint Miss Standard Deviation (Ft)");
    let series1 = vec![
        Series::new(results.tp.clone(), results.xmgl.clone())
            .with_label("Simulation")
            .with_color(plotters::prelude::BLUE),
        Series::new(results.tp.clone(), results.xmgl_th.clone())
            .with_label("Theory")
            .with_color(plotters::prelude::RED),
    ];
    line_plot(&plot_file1, &config1, &series1).ok();

    let plot_file2 = format!("{}/c4l7_range_indep.png", output_dir);
    let config2 = PlotConfig::new("Range Independent Noise Miss Standard Deviation")
        .with_labels("Flight Time (Sec)", "Miss Standard Deviation (Ft)");
    let series2 = vec![
        Series::new(results.tp.clone(), results.xmfn.clone())
            .with_label("Simulation")
            .with_color(plotters::prelude::BLUE),
        Series::new(results.tp.clone(), results.xmfn_th.clone())
            .with_label("Theory")
            .with_color(plotters::prelude::RED),
    ];
    line_plot(&plot_file2, &config2, &series2).ok();

    let plot_file3 = format!("{}/c4l7_semiactive.png", output_dir);
    let config3 = PlotConfig::new("Semiactive Range Dependent Noise Miss")
        .with_labels("Flight Time (Sec)", "Miss Standard Deviation (Ft)");
    let series3 = vec![
        Series::new(results.tp.clone(), results.xmrn.clone())
            .with_label("Simulation")
            .with_color(plotters::prelude::BLUE),
        Series::new(results.tp.clone(), results.xmrn_th.clone())
            .with_label("Theory")
            .with_color(plotters::prelude::RED),
    ];
    line_plot(&plot_file3, &config3, &series3).ok();

    let plot_file4 = format!("{}/c4l7_active.png", output_dir);
    let config4 = PlotConfig::new("Active Range Dependent Noise Miss")
        .with_labels("Flight Time (Sec)", "Miss Standard Deviation (Ft)");
    let series4 = vec![
        Series::new(results.tp.clone(), results.xmrna.clone())
            .with_label("Simulation")
            .with_color(plotters::prelude::BLUE),
        Series::new(results.tp.clone(), results.xmrna_th.clone())
            .with_label("Theory")
            .with_color(plotters::prelude::RED),
    ];
    line_plot(&plot_file4, &config4, &series4).ok();

    println!("C4L7: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c4l7_runs() {
        let results = run();
        assert!(results.tp.len() > 0);
        assert_eq!(results.tp.len(), results.xmgl.len());
    }
}
