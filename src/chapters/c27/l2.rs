//! Chapter 27, Lesson 2: Target Maneuver Miss Analysis with Fading Memory Filter
//!
//! Uses fading memory filter for LOS rate estimation with digital filter
//! for command smoothing.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub tp: Vec<f64>,
    pub xmnt: Vec<f64>,
}

/// Run the C27L2 simulation
pub fn run() -> Results {
    let xnt: f64 = 96.6;
    let xnp: f64 = 3.0;
    let tf: f64 = 10.0;
    let ts: f64 = 0.1;
    let beta: f64 = 0.8;
    let tap: f64 = 0.2;
    let vc: f64 = 4000.0;
    let ts2: f64 = 0.02;
    let h: f64 = 0.001;

    let gfilter = 1.0 - beta * beta;
    let hfilter = (1.0 - beta).powi(2);

    let mut tp = 0.00001;
    let mut s: f64 = 0.0;
    let mut s2: f64 = 0.0;

    let mut x1: f64 = 0.0;
    let mut x2: f64 = 0.0;
    let mut x3: f64 = 1.0;
    let mut x4: f64 = 0.0;
    let mut x5: f64 = 0.0;

    let mut y1old: f64 = 0.0;
    let mut y2old: f64 = 0.0;
    let mut y3old: f64 = 0.0;
    let mut y4old: f64 = 0.0;
    let mut y5old: f64 = 0.0;
    let mut y6old: f64 = 0.0;
    let mut y7old: f64 = 0.0;
    let mut y8old: f64 = 0.0;
    let mut y9old: f64 = 0.0;
    let mut y10old: f64 = 0.0;
    let mut y11old: f64 = 0.0;

    let mut y1new: f64;
    let mut y2new: f64;
    let mut y3new: f64;
    let mut y4new: f64 = 0.0;
    let mut y5new: f64;
    let mut y6new: f64;
    let mut y7new: f64 = 0.0;
    let mut y8new: f64;
    let mut y9new: f64;
    let mut y10new: f64;
    let mut y11new: f64;

    let mut array_tp = Vec::new();
    let mut array_xmnt = Vec::new();

    // MATLAB initializes count=1 then increments to 2 before storing,
    // leaving index 1 as implicit zeros. We match this by starting with zeros.
    array_tp.push(0.0);
    array_xmnt.push(0.0);

    while tp <= tf - 1e-5 {
        s += h;
        s2 += h;

        let x1old = x1;
        let x2old = x2;
        let x3old = x3;
        let x4old = x4;
        let x5old = x5;

        // First derivative evaluation (FLAG=0 initially)
        let tgo = tp;
        let x1d = x2;
        let x2d = x3 + y4new / (vc * tgo);
        let x3d = y4new / (vc * tgo * tgo);
        let x4d = -x2 - x4 / tap;
        let x5d = x4 / tap;

        // Euler step (FLAG=1)
        x1 += h * x1d;
        x2 += h * x2d;
        x3 += h * x3d;
        x4 += h * x4d;
        x5 += h * x5d;
        tp += h;

        // Second derivative for RK2
        let tgo = tp;
        let x1d = x2;
        let x2d = x3 + y4new / (vc * tgo);
        let x3d = y4new / (vc * tgo * tgo);
        let x4d = -x2 - x4 / tap;
        let x5d = x4 / tap;

        // RK2 averaging
        x1 = (x1old + x1) / 2.0 + 0.5 * h * x1d;
        x2 = (x2old + x2) / 2.0 + 0.5 * h * x2d;
        x3 = (x3old + x3) / 2.0 + 0.5 * h * x3d;
        x4 = (x4old + x4) / 2.0 + 0.5 * h * x4d;
        x5 = (x5old + x5) / 2.0 + 0.5 * h * x5d;

        // TS sampling period for fading memory filter
        if s >= ts - 0.0001 {
            s = 0.0;
            y1new = x5;
            let temp1 = (y1new - y1old) * xnp * vc;
            let temp2 = hfilter * (y2old + temp1) / ts + gfilter * y3old;
            y2new = temp1 + y2old + ts * (y3old - temp2);
            y3new = y3old - temp2;
            y7new = temp2 + y7old;
            y1old = y1new;
            y2old = y2new;
            y3old = y3new;
            y7old = y7new;
        }

        // TS2 sampling period for digital filter
        if s2 >= ts2 - 0.0001 {
            s2 = 0.0;
            y6new = y7new;
            y11new = 0.2 * (y6new - y6old);
            y10new = y11old + y11new;
            y9new = y10old + y11new;
            y8new = y9old + y11new;
            y5new = y8old + y11new;
            y4new = y4old + y5old + y11new;
            y4old = y4new;
            y6old = y6new;
            y5old = y5new;
            y8old = y8new;
            y9old = y9new;
            y10old = y10new;
            y11old = y11new;
            let xmnt = xnt * x1;
            array_tp.push(tp);
            array_xmnt.push(xmnt);
        }
    }

    Results {
        tp: array_tp,
        xmnt: array_xmnt,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c27l2_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.tp.clone(),
        results.xmnt.clone(),
    ])?;

    let plot_file = format!("{}/c27l2_miss.png", output_dir);
    let config = PlotConfig::new("Target Maneuver Miss")
        .with_labels("Flight Time (Sec)", "Target Maneuver Miss (Ft)");

    let series = vec![
        Series::new(results.tp.clone(), results.xmnt.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C27L2: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c27l2_runs() {
        let results = run();
        assert!(!results.tp.is_empty());
    }
}
