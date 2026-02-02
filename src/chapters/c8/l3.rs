//! Chapter 8, Lesson 3: Optimal Proportional Navigation with TGO Error
//!
//! Adjoint analysis for OPN with time-to-go estimation error.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub xmnt: Vec<f64>,
    pub xmhe: Vec<f64>,
}

/// Run the C8L3 simulation
pub fn run() -> Results {
    let xnt: f64 = 96.6;
    let xnp: f64 = 4.0;
    let tau: f64 = 1.0;
    let tf: f64 = 10.0;
    let vm: f64 = 3000.0;
    let hedeg: f64 = -20.0;
    let apn: i32 = 2;  // 0 = PN, 1 = APN, 2 = OPN
    let bias: f64 = 0.1;
    let sf: f64 = 1.0;
    let h: f64 = 0.01;
    let he = hedeg / 57.3;

    let mut tp: f64 = 0.00001;
    let mut s: f64 = 0.0;

    let mut x1: f64 = 0.0;
    let mut x2: f64 = 0.0;
    let mut x3: f64 = 1.0;
    let mut x4: f64 = 0.0;

    let mut array_tp = Vec::new();
    let mut array_xmnt = Vec::new();
    let mut array_xmhe = Vec::new();

    while tp <= (tf - 1e-5) {
        let x1_old = x1;
        let x2_old = x2;
        let x3_old = x3;
        let x4_old = x4;

        // First derivative evaluation
        let tgo = tp + 0.00001;
        let (c1, c2, c3, c4) = match apn {
            0 => (xnp / (tgo * tgo), xnp / tgo, 0.0, 0.0),
            1 => (xnp / (tgo * tgo), xnp / tgo, 0.5 * xnp, 0.0),
            _ => {
                let mut tgoh = sf * tgo + bias;
                if tgoh < 0.0 {
                    tgoh = 0.0001;
                }
                let x = tgoh / tau;
                let top = 6.0 * x * x * ((-x).exp() - 1.0 + x);
                let bot1 = 2.0 * x * x * x + 3.0 + 6.0 * x - 6.0 * x * x;
                let bot2 = -12.0 * x * (-x).exp() - 3.0 * (-2.0 * x).exp();
                let xnpp = top / (0.0001 + bot1 + bot2);
                (
                    xnpp / (tgoh * tgoh),
                    xnpp / tgoh,
                    0.5 * xnpp,
                    -xnpp * ((-x).exp() + x - 1.0) / (x * x),
                )
            }
        };

        let x1d = x2 + c3 * x4 / tau;
        let x2d = x3 + c2 * x4 / tau;
        let x3d = c1 * x4 / tau;
        let x4d = -x4 / tau - x2 + c4 * x4 / tau;

        // Euler step
        x1 += h * x1d;
        x2 += h * x2d;
        x3 += h * x3d;
        x4 += h * x4d;
        tp += h;

        // Second derivative
        let tgo = tp + 0.00001;
        let (c1, c2, c3, c4) = match apn {
            0 => (xnp / (tgo * tgo), xnp / tgo, 0.0, 0.0),
            1 => (xnp / (tgo * tgo), xnp / tgo, 0.5 * xnp, 0.0),
            _ => {
                let mut tgoh = sf * tgo + bias;
                if tgoh < 0.0 {
                    tgoh = 0.0001;
                }
                let x = tgoh / tau;
                let top = 6.0 * x * x * ((-x).exp() - 1.0 + x);
                let bot1 = 2.0 * x * x * x + 3.0 + 6.0 * x - 6.0 * x * x;
                let bot2 = -12.0 * x * (-x).exp() - 3.0 * (-2.0 * x).exp();
                let xnpp = top / (0.0001 + bot1 + bot2);
                (
                    xnpp / (tgoh * tgoh),
                    xnpp / tgoh,
                    0.5 * xnpp,
                    -xnpp * ((-x).exp() + x - 1.0) / (x * x),
                )
            }
        };

        let x1d = x2 + c3 * x4 / tau;
        let x2d = x3 + c2 * x4 / tau;
        let x3d = c1 * x4 / tau;
        let x4d = -x4 / tau - x2 + c4 * x4 / tau;

        // RK2 averaging
        x1 = (x1_old + x1) / 2.0 + 0.5 * h * x1d;
        x2 = (x2_old + x2) / 2.0 + 0.5 * h * x2d;
        x3 = (x3_old + x3) / 2.0 + 0.5 * h * x3d;
        x4 = (x4_old + x4) / 2.0 + 0.5 * h * x4d;

        s += h;
        if s >= 0.0999 {
            s = 0.0;
            array_tp.push(tp);
            array_xmnt.push(xnt * x1);
            array_xmhe.push(-vm * he * x2);
        }
    }

    Results {
        time: array_tp,
        xmnt: array_xmnt,
        xmhe: array_xmhe,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c8l3_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.xmnt.clone(),
        results.xmhe.clone(),
    ])?;

    let plot_file = format!("{}/c8l3_plot.png", output_dir);
    let config = PlotConfig::new("Target Maneuver Miss (OPN with TGO Error)")
        .with_labels("Flight Time (Sec)", "Target Maneuver Miss (Ft)");
    let series = vec![
        Series::new(results.time.clone(), results.xmnt.clone())
            .with_color(plotters::prelude::BLUE),
    ];
    line_plot(&plot_file, &config, &series).ok();

    println!("C8L3: Simulation finished");
    println!("  Data saved to: {}", data_file);
    println!("  Plot saved to: {}", plot_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c8l3_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
