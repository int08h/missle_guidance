//! Chapter 7, Lesson 3: Digital Fading Memory Filter Adjoint
//!
//! Adjoint analysis for digital fading memory filter.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub xmnt: Vec<f64>,
    pub xmnoise: Vec<f64>,
}

/// Run the C7L3 simulation
pub fn run() -> Results {
    let xnt: f64 = 96.6;
    let xnp: f64 = 3.0;
    let tf: f64 = 10.0;
    let ts: f64 = 0.1;
    let beta: f64 = 0.8;
    let signoise: f64 = 0.001;
    let vc: f64 = 4000.0;
    let h: f64 = 0.01;

    let gfilter = 1.0 - beta * beta;
    let hfilter = (1.0 - beta).powi(2);

    let mut tp: f64 = 0.00001;
    let mut s: f64 = 0.0;

    let mut x1: f64 = 0.0;
    let mut x2: f64 = 0.0;
    let mut x3: f64 = 1.0;
    let mut x5: f64 = 0.0;

    let mut y1_old: f64 = 0.0;
    let mut y2_old: f64 = 0.0;
    let mut y3_old: f64 = 0.0;
    let mut y4_old: f64 = 0.0;
    let mut y5_old: f64 = 0.0;

    let mut array_tp = Vec::new();
    let mut array_xmnt = Vec::new();
    let mut array_xmnoise = Vec::new();

    // First entry at count=1 is skipped (MATLAB starts count at 1 then increments)

    while tp <= (tf - 1e-5) {
        let x1_old = x1;
        let x2_old = x2;
        let x3_old = x3;
        let x5_old = x5;

        // First derivative evaluation
        let tgo = tp + 0.00001;
        let x1d = x2;
        let x2d = x3 + y4_old / (vc * tgo);
        let x3d = y4_old / (vc * tgo * tgo);
        let x5d = -x2;

        // Euler step
        x1 += h * x1d;
        x2 += h * x2d;
        x3 += h * x3d;
        x5 += h * x5d;
        tp += h;

        // Second derivative
        let tgo = tp + 0.00001;
        let x1d = x2;
        let x2d = x3 + y4_old / (vc * tgo);
        let x3d = y4_old / (vc * tgo * tgo);
        let x5d = -x2;

        // RK2 averaging
        x1 = (x1_old + x1) / 2.0 + 0.5 * h * x1d;
        x2 = (x2_old + x2) / 2.0 + 0.5 * h * x2d;
        x3 = (x3_old + x3) / 2.0 + 0.5 * h * x3d;
        x5 = (x5_old + x5) / 2.0 + 0.5 * h * x5d;

        s += h;
        if s >= (ts - 0.0001) {
            s = 0.0;
            let temp1 = (x5 - y1_old) * xnp * vc;
            let temp2 = hfilter * (y2_old + temp1) / ts + gfilter * y3_old;
            let y1_new = x5;
            let y2_new = temp1 + y2_old + ts * (y3_old - temp2);
            let y3_new = y3_old - temp2;
            let y4_new = y4_old + temp2;
            let y5_new = y5_old + temp2 * temp2;

            y1_old = y1_new;
            y2_old = y2_new;
            y3_old = y3_new;
            y4_old = y4_new;
            y5_old = y5_new;

            let xmnoise = signoise * y5_new.sqrt();
            let xmnt = xnt * x1;

            array_tp.push(tp);
            array_xmnt.push(xmnt);
            array_xmnoise.push(xmnoise);
        }
    }

    Results {
        time: array_tp,
        xmnt: array_xmnt,
        xmnoise: array_xmnoise,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c7l3_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.xmnt.clone(),
        results.xmnoise.clone(),
    ])?;

    let plot_file = format!("{}/c7l3_maneuver.png", output_dir);
    let config = PlotConfig::new("Target Maneuver Miss")
        .with_labels("Flight Time (Sec)", "Target Maneuver Miss (Ft)");
    let series = vec![
        Series::new(results.time.clone(), results.xmnt.clone())
            .with_color(plotters::prelude::BLUE),
    ];
    line_plot(&plot_file, &config, &series).ok();

    let plot_file2 = format!("{}/c7l3_noise.png", output_dir);
    let config2 = PlotConfig::new("Noise Miss")
        .with_labels("Flight Time (Sec)", "Noise Miss (Ft)");
    let series2 = vec![
        Series::new(results.time.clone(), results.xmnoise.clone())
            .with_color(plotters::prelude::BLUE),
    ];
    line_plot(&plot_file2, &config2, &series2).ok();

    println!("C7L3: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c7l3_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
