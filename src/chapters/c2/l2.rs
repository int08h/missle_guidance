//! Chapter 2, Lesson 2: Linearized Engagement Model
//!
//! Simulates proportional navigation guidance using a linearized model.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

/// Simulation results
pub struct Results {
    pub time: Vec<f64>,
    pub y: Vec<f64>,
    pub yd: Vec<f64>,
    pub xnc_g: Vec<f64>,
}

/// Run the C2L2 simulation
pub fn run() -> Results {
    let vc = 4000.0;     // Closing velocity (ft/s)
    let xnt = 0.0;       // Target acceleration
    let vm = 3000.0;     // Missile velocity (ft/s)
    let he_deg = -20.0;  // Heading error (degrees)
    let tf = 10.0;       // Flight time
    let xnp = 4.0;       // Navigation ratio

    let mut y = 0.0;
    let mut yd = -vm * he_deg / 57.3;
    let mut t = 0.0;
    let h = 0.01;
    let mut s = 0.0;

    let mut xnc: f64;

    let mut array_t = Vec::new();
    let mut array_y = Vec::new();
    let mut array_yd = Vec::new();
    let mut array_xnc_g = Vec::new();

    while t <= tf - 1e-5 {
        // Save old values
        let yold = y;
        let ydold = yd;

        // First derivative evaluation
        let tgo = tf - t + 0.00001;
        let xlamd = (y + yd * tgo) / (vc * tgo * tgo);
        xnc = xnp * vc * xlamd;
        let ydd = xnt - xnc;

        // Euler step
        y += h * yd;
        yd += h * ydd;
        t += h;

        // Second evaluation and RK2 averaging
        let tgo = tf - t + 0.00001;
        let xlamd = (y + yd * tgo) / (vc * tgo * tgo);
        xnc = xnp * vc * xlamd;
        let ydd = xnt - xnc;

        y = 0.5 * (yold + y + h * yd);
        yd = 0.5 * (ydold + yd + h * ydd);

        s += h;

        // Store data at sampling interval
        if s >= 0.0999 {
            s = 0.0;
            array_t.push(t);
            array_y.push(y);
            array_yd.push(yd);
            array_xnc_g.push(xnc / 32.2);
        }
    }

    Results {
        time: array_t,
        y: array_y,
        yd: array_yd,
        xnc_g: array_xnc_g,
    }
}

/// Run and save results to file
pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    // Save data file
    let data_file = format!("{}/c2l2_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.y.clone(),
        results.yd.clone(),
        results.xnc_g.clone(),
    ])?;

    // Create acceleration plot
    let plot_file = format!("{}/c2l2_plot.png", output_dir);
    let config = PlotConfig::new("Linearized PN Guidance")
        .with_labels("Time (Sec)", "Missile Acceleration (G)");

    let series = vec![
        Series::new(results.time.clone(), results.xnc_g.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C2L2: Simulation finished");
    println!("  Data saved to: {}", data_file);
    println!("  Plot saved to: {}", plot_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c2l2_runs() {
        let results = run();
        assert!(!results.time.is_empty());
        assert_eq!(results.time.len(), results.xnc_g.len());
    }
}
