//! Chapter 35, Lesson 3: Adjoint with First-Order Lag
//!
//! Computes optimal guidance gains using adjoint method
//! with first-order lag autopilot dynamics.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub c4: Vec<f64>,
    pub c4th: Vec<f64>,
    pub xnp: Vec<f64>,
    pub xnpth: Vec<f64>,
}

/// Run the C35L3 simulation
pub fn run() -> Results {
    let tau: f64 = 1.0;
    let gam: f64 = 0.00001;

    let mut t: f64 = 0.0;
    let h: f64 = 0.01;
    let mut s: f64 = 0.0;

    // State variables
    let mut x1: f64 = 0.0;
    let mut x2: f64 = 0.0;
    let mut x3: f64 = 0.0;

    let mut array_t = Vec::new();
    let mut array_c4 = Vec::new();
    let mut array_c4th = Vec::new();
    let mut array_xnp = Vec::new();
    let mut array_xnpth = Vec::new();

    while t < 10.0 - 0.0001 {
        s += h;
        let x1old = x1;
        let x2old = x2;
        let x3old = x3;

        // First derivative evaluation
        let x1d = (t - x1) / tau;
        let x2d = x1 * x1;
        let _d = x2 + gam;
        let x3d = -x3 / tau + t;

        // Euler step
        x1 += h * x1d;
        x2 += h * x2d;
        x3 += h * x3d;
        t += h;

        // Second derivative for RK2
        let x1d = (t - x1) / tau;
        let x2d = x1 * x1;
        let d = x2 + gam;
        let x3d = -x3 / tau + t;

        // RK2 averaging
        x1 = 0.5 * (x1old + x1 + h * x1d);
        x2 = 0.5 * (x2old + x2 + h * x2d);
        x3 = 0.5 * (x3old + x3 + h * x3d);

        if s >= 0.09999 {
            s = 0.0;
            let pz = x1 / d;
            let xnp = pz * t * t;
            let _c1 = pz;
            let _c2 = pz * t;
            let _c3 = 0.5 * pz * t * t;
            let c4 = -x3 * pz;

            // Theoretical values
            let xs = t / tau;
            let top = 6.0 * xs * xs * ((-xs).exp() - 1.0 + xs);
            let bot1 = 2.0 * xs.powi(3) + 3.0 + 6.0 * xs - 6.0 * xs * xs;
            let bot2 = -12.0 * xs * (-xs).exp() - 3.0 * (-2.0 * xs).exp();
            let xnpth = top / (bot1 + bot2 + 0.0001);
            let c4th = -xnpth * ((-xs).exp() + xs - 1.0) / (xs * xs);

            array_t.push(t);
            array_c4.push(c4);
            array_c4th.push(c4th);
            array_xnp.push(xnp);
            array_xnpth.push(xnpth);
        }
    }

    Results {
        time: array_t,
        c4: array_c4,
        c4th: array_c4th,
        xnp: array_xnp,
        xnpth: array_xnpth,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c35l3_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.c4.clone(),
        results.c4th.clone(),
        results.xnp.clone(),
        results.xnpth.clone(),
    ])?;

    let plot_file = format!("{}/c35l3_np.png", output_dir);
    let config = PlotConfig::new("Navigation Ratio NP")
        .with_labels("Time (s)", "NP");

    let series = vec![
        Series::new(results.time.clone(), results.xnp.clone())
            .with_color(plotters::prelude::BLUE)
            .with_label("NP"),
        Series::new(results.time.clone(), results.xnpth.clone())
            .with_color(plotters::prelude::RED)
            .with_label("NP Theory"),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C35L3: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c35l3_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
