//! Chapter 4, Lesson 8: Error Budget Analysis
//!
//! Computes RMS miss for semiactive and active error budgets combining
//! glint, range-independent noise, range-dependent noise, and target maneuver.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

/// Simulation results
pub struct Results {
    pub tp: Vec<f64>,
    pub xmgl: Vec<f64>,
    pub xmfn: Vec<f64>,
    pub xmrn: Vec<f64>,
    pub xmrna: Vec<f64>,
    pub rmssa: Vec<f64>,
    pub rmsa: Vec<f64>,
}

/// Run the C4L8 simulation
pub fn run() -> Results {
    let xnt: f64 = 193.2;
    let xnp: f64 = 3.0;
    let tau: f64 = 0.5;
    let tf: f64 = 10.0;
    let vc: f64 = 4000.0;
    let ts: f64 = 0.1;
    let ra: f64 = 30000.0;

    // Noise parameters
    let siggl: f64 = 2.0;
    let sigfn: f64 = 0.001;
    let sigrn: f64 = 0.01;
    let sigrna: f64 = 0.01;

    let phigl: f64 = siggl * siggl / ts;
    let phifn: f64 = sigfn * sigfn / ts;
    let phirn: f64 = sigrn * sigrn / ts;
    let phirna: f64 = sigrna * sigrna / ts;

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
    let mut x9 = 0.0;

    let mut array_tp = Vec::new();
    let mut array_xmgl = Vec::new();
    let mut array_xmfn = Vec::new();
    let mut array_xmrn = Vec::new();
    let mut array_xmrna = Vec::new();
    let mut array_rmssa = Vec::new();
    let mut array_rmsa = Vec::new();

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
        let x9_old = x9;

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
        let x9d = x1 * x1;

        x1 += h * x1d;
        x2 += h * x2d;
        x3 += h * x3d;
        x4 += h * x4d;
        x5 += h * x5d;
        x6 += h * x6d;
        x7 += h * x7d;
        x8 += h * x8d;
        x9 += h * x9d;
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
        let x9d = x1 * x1;

        x1 = (x1_old + x1) / 2.0 + 0.5 * h * x1d;
        x2 = (x2_old + x2) / 2.0 + 0.5 * h * x2d;
        x3 = (x3_old + x3) / 2.0 + 0.5 * h * x3d;
        x4 = (x4_old + x4) / 2.0 + 0.5 * h * x4d;
        x5 = (x5_old + x5) / 2.0 + 0.5 * h * x5d;
        x6 = (x6_old + x6) / 2.0 + 0.5 * h * x6d;
        x7 = (x7_old + x7) / 2.0 + 0.5 * h * x7d;
        x8 = (x8_old + x8) / 2.0 + 0.5 * h * x8d;
        x9 = (x9_old + x9) / 2.0 + 0.5 * h * x9d;

        s += h;
        if s >= 0.0999 {
            s = 0.0;
            let xmgl = (phigl * x5).sqrt();
            let xmfn = (phifn * x6).sqrt();
            let xmrn = (phirn * x7).sqrt();
            let xmrna = (phirna * x8).sqrt();
            let xmudnt = xnt * (x9 / tgo).sqrt();

            // RMS miss for semiactive (uses xmrn)
            let rmssa = (xmgl.powi(2) + xmfn.powi(2) + xmrn.powi(2) + xmudnt.powi(2)).sqrt();
            // RMS miss for active (uses xmrna)
            let rmsa = (xmgl.powi(2) + xmfn.powi(2) + xmrna.powi(2) + xmudnt.powi(2)).sqrt();

            array_tp.push(tp);
            array_xmgl.push(xmgl);
            array_xmfn.push(xmfn);
            array_xmrn.push(xmrn);
            array_xmrna.push(xmrna);
            array_rmssa.push(rmssa);
            array_rmsa.push(rmsa);
        }
    }

    Results {
        tp: array_tp,
        xmgl: array_xmgl,
        xmfn: array_xmfn,
        xmrn: array_xmrn,
        xmrna: array_xmrna,
        rmssa: array_rmssa,
        rmsa: array_rmsa,
    }
}

/// Run and save results to file
pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    // Save data file
    let data_file = format!("{}/c4l8_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.tp.clone(),
        results.xmgl.clone(),
        results.xmfn.clone(),
        results.xmrn.clone(),
        results.xmrna.clone(),
        results.rmssa.clone(),
        results.rmsa.clone(),
    ])?;

    // Create plots
    let plot_file1 = format!("{}/c4l8_semiactive.png", output_dir);
    let config1 = PlotConfig::new("RMS Miss For Semiactive Error Budget")
        .with_labels("Flight Time (Sec)", "RMS Miss (Ft)");
    let series1 = vec![
        Series::new(results.tp.clone(), results.xmgl.clone())
            .with_label("Glint")
            .with_color(plotters::prelude::BLUE),
        Series::new(results.tp.clone(), results.xmfn.clone())
            .with_label("FN")
            .with_color(plotters::prelude::RED),
        Series::new(results.tp.clone(), results.xmrn.clone())
            .with_label("RN")
            .with_color(plotters::prelude::GREEN),
        Series::new(results.tp.clone(), results.rmssa.clone())
            .with_label("Total")
            .with_color(plotters::prelude::BLACK),
    ];
    line_plot(&plot_file1, &config1, &series1).ok();

    let plot_file2 = format!("{}/c4l8_active.png", output_dir);
    let config2 = PlotConfig::new("RMS Miss For Active Error Budget")
        .with_labels("Flight Time (Sec)", "RMS Miss (Ft)");
    let series2 = vec![
        Series::new(results.tp.clone(), results.xmgl.clone())
            .with_label("Glint")
            .with_color(plotters::prelude::BLUE),
        Series::new(results.tp.clone(), results.xmfn.clone())
            .with_label("FN")
            .with_color(plotters::prelude::RED),
        Series::new(results.tp.clone(), results.xmrna.clone())
            .with_label("RNA")
            .with_color(plotters::prelude::GREEN),
        Series::new(results.tp.clone(), results.rmsa.clone())
            .with_label("Total")
            .with_color(plotters::prelude::BLACK),
    ];
    line_plot(&plot_file2, &config2, &series2).ok();

    println!("C4L8: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c4l8_runs() {
        let results = run();
        assert!(results.tp.len() > 0);
        assert_eq!(results.tp.len(), results.xmgl.len());
    }
}
