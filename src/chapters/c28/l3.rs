//! Chapter 28, Lesson 3: Adjoint Model with Shaping Filter
//!
//! Uses adjoint method with shaping filter approach to compute
//! miss distance standard deviation.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub tp: Vec<f64>,
    pub xmudnt: Vec<f64>,
}

/// Run the C28L3 simulation
pub fn run() -> Results {
    let xnt: f64 = 161.0;
    let xnp: f64 = 3.0;
    let tau: f64 = 0.2;
    let tf: f64 = 10.0;
    let h: f64 = 0.01;

    let mut tp = 0.00001;
    let mut s: f64 = 0.0;

    let mut x1: f64 = 0.0;
    let mut x2: f64 = 0.0;
    let mut x3: f64 = 1.0;
    let mut x4: f64 = 0.0;
    let mut x5: f64 = 0.0;

    let mut array_tp = Vec::new();
    let mut array_xmudnt = Vec::new();

    while tp <= tf - 1e-5 {
        s += h;

        let x1old = x1;
        let x2old = x2;
        let x3old = x3;
        let x4old = x4;
        let x5old = x5;

        // First derivative evaluation
        let x1d = x2;
        let x2d = x3;
        let y1 = (x4 - x2) / tau;
        let tgo = tp + 0.00001;
        let x3d = xnp * y1 / tgo;
        let x4d = -y1;
        let x5d = x1 * x1;

        // Euler step
        x1 += h * x1d;
        x2 += h * x2d;
        x3 += h * x3d;
        x4 += h * x4d;
        x5 += h * x5d;
        tp += h;

        // Second derivative for RK2
        let x1d = x2;
        let x2d = x3;
        let y1 = (x4 - x2) / tau;
        let tgo = tp + 0.00001;
        let x3d = xnp * y1 / tgo;
        let x4d = -y1;
        let x5d = x1 * x1;

        // RK2 averaging
        x1 = (x1old + x1) / 2.0 + 0.5 * h * x1d;
        x2 = (x2old + x2) / 2.0 + 0.5 * h * x2d;
        x3 = (x3old + x3) / 2.0 + 0.5 * h * x3d;
        x4 = (x4old + x4) / 2.0 + 0.5 * h * x4d;
        x5 = (x5old + x5) / 2.0 + 0.5 * h * x5d;

        // Additional s increment per MATLAB line 50
        s += h;

        if s >= 0.000999 {
            s = 0.0;
            let tgo = tp + 0.00001;
            let xmudnt = xnt * (x5 / tgo).sqrt();
            array_tp.push(tp);
            array_xmudnt.push(xmudnt);
        }
    }

    Results {
        tp: array_tp,
        xmudnt: array_xmudnt,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c28l3_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.tp.clone(),
        results.xmudnt.clone(),
    ])?;

    let plot_file = format!("{}/c28l3_miss.png", output_dir);
    let config = PlotConfig::new("Adjoint Model Using Shaping Filter Approach")
        .with_labels("Flight Time (S)", "Miss Dist Standard Deviation (Ft)");

    let series = vec![
        Series::new(results.tp.clone(), results.xmudnt.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C28L3: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c28l3_runs() {
        let results = run();
        assert!(!results.tp.is_empty());
    }
}
