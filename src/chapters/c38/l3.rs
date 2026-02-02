//! Chapter 38, Lesson 3: Adjoint Model with Trapezoidal Weave
//!
//! Computes miss distance using adjoint method with trapezoidal weave
//! target maneuver model.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

// Use MATLAB's PI value (3.1416) instead of std::f64::consts::PI
// to ensure numerical compatibility with the original code.
const PI: f64 = 3.1416;

pub struct Results {
    pub time: Vec<f64>,
    pub xmnt: Vec<f64>,
    pub xmudnt: Vec<f64>,
}

/// Run the C38L3 simulation
pub fn run() -> Results {
    let xnt: f64 = 161.0;
    let xnp: f64 = 3.0;
    let tau: f64 = 0.5;
    let tf: f64 = 10.0;
    let _vm: f64 = 2000.0;
    let _hedeg: f64 = 0.0;
    let pz: f64 = 3.0;
    let xl = pz / 2.0;
    let tr: f64 = 1.0;
    let alf = PI * tr / (2.0 * xl);
    let w = 2.0 * PI / pz;
    let vc: f64 = 3000.0;

    let t: f64 = 0.0;
    let mut s: f64 = 0.0;
    let mut tp = t + 0.00001;

    // State variables
    let mut x1: f64 = 0.0;
    let mut x2: f64 = 0.0;
    let mut x3: f64 = 1.0;
    let mut x4: f64 = 0.0;
    let mut x5: f64 = 0.0;
    let mut x6: f64 = 0.0;
    let mut x7: f64 = 0.0;
    let mut x8: f64 = 0.0;

    let h: f64 = 0.01;
    let _he = _hedeg / 57.3;

    let mut array_tp = Vec::new();
    let mut array_xmnt = Vec::new();
    let mut array_xmudnt = Vec::new();

    while tp <= tf - 1e-5 {
        let x1old = x1;
        let x2old = x2;
        let x3old = x3;
        let x4old = x4;
        let x5old = x5;
        let x6old = x6;
        let x7old = x7;
        let x8old = x8;

        // First derivative evaluation
        let tgo = tp + 0.00001;
        let x1d = x5;
        let x2d = x3;
        let y1 = x4 - xnp * vc * x2;
        let x3d = y1 / (vc * tgo * tau);
        let x4d = -y1 / tau;
        let x5d = x2 - w * w * x1;
        let x6d = x7;
        let x7d = x2 - 9.0 * w * w * x6;
        let y2 = x1 * alf.sin() + x6 * (3.0 * alf).sin() / 9.0;
        let pz_val = 4.0 * w * y2 / (PI * alf);
        let x8d = pz_val * pz_val;

        // Euler step
        x1 += h * x1d;
        x2 += h * x2d;
        x3 += h * x3d;
        x4 += h * x4d;
        x5 += h * x5d;
        x6 += h * x6d;
        x7 += h * x7d;
        x8 += h * x8d;
        tp += h;

        // Second derivative for RK2
        let tgo = tp + 0.00001;
        let x1d = x5;
        let x2d = x3;
        let y1 = x4 - xnp * vc * x2;
        let x3d = y1 / (vc * tgo * tau);
        let x4d = -y1 / tau;
        let x5d = x2 - w * w * x1;
        let x6d = x7;
        let x7d = x2 - 9.0 * w * w * x6;
        let y2 = x1 * alf.sin() + x6 * (3.0 * alf).sin() / 9.0;
        let pz_val = 4.0 * w * y2 / (PI * alf);
        let x8d = pz_val * pz_val;

        // RK2 averaging
        x1 = (x1old + x1) / 2.0 + 0.5 * h * x1d;
        x2 = (x2old + x2) / 2.0 + 0.5 * h * x2d;
        x3 = (x3old + x3) / 2.0 + 0.5 * h * x3d;
        x4 = (x4old + x4) / 2.0 + 0.5 * h * x4d;
        x5 = (x5old + x5) / 2.0 + 0.5 * h * x5d;
        x6 = (x6old + x6) / 2.0 + 0.5 * h * x6d;
        x7 = (x7old + x7) / 2.0 + 0.5 * h * x7d;
        x8 = (x8old + x8) / 2.0 + 0.5 * h * x8d;

        s += h;
        if s >= 0.0999 {
            s = 0.0;
            let xmnt = 4.0 * w * xnt * y2 / (PI * alf);
            let xmudnt = xnt * (x8 / tgo).sqrt();

            array_tp.push(tp);
            array_xmnt.push(xmnt);
            array_xmudnt.push(xmudnt);
        }
    }

    Results {
        time: array_tp,
        xmnt: array_xmnt,
        xmudnt: array_xmudnt,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c38l3_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.xmnt.clone(),
        results.xmudnt.clone(),
    ])?;

    let plot_file = format!("{}/c38l3_miss.png", output_dir);
    let config = PlotConfig::new("Target Maneuver Miss")
        .with_labels("Flight Time (Sec)", "Miss (Ft)");

    let series = vec![
        Series::new(results.time.clone(), results.xmnt.clone())
            .with_color(plotters::prelude::BLUE)
            .with_label("Deterministic"),
        Series::new(results.time.clone(), results.xmudnt.clone())
            .with_color(plotters::prelude::RED)
            .with_label("Random"),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C38L3: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c38l3_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
