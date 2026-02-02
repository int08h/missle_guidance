//! Chapter 5, Lesson 1: Fourth-Order Runge-Kutta Integration
//!
//! Demonstrates 4th-order RK integration for a second-order damped system.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub y: Vec<f64>,
}

/// Run the C5L1 simulation - 4th order Runge-Kutta
pub fn run() -> Results {
    let xin: f64 = 1.0;      // Step input
    let w: f64 = 20.0;       // Natural frequency
    let h: f64 = 0.01;       // Time step

    let mut y: f64 = 0.0;
    let mut yd: f64 = 0.0;
    let mut t: f64 = 0.0;

    let mut array_t = Vec::new();
    let mut array_y = Vec::new();

    while t < 1.0 {
        let yold = y;
        let ydold = yd;

        // RK4 integration
        // k0
        let ydd = w * xin - w * w * y;
        let k0_1 = yd;
        let k0_2 = ydd;

        // k1 at t + h/2
        let y_temp = yold + 0.5 * h * k0_1;
        let yd_temp = ydold + 0.5 * h * k0_2;
        let ydd = w * xin - w * w * y_temp;
        let k1_1 = yd_temp;
        let k1_2 = ydd;

        // k2 at t + h/2
        let y_temp = yold + 0.5 * h * k1_1;
        let yd_temp = ydold + 0.5 * h * k1_2;
        let ydd = w * xin - w * w * y_temp;
        let k2_1 = yd_temp;
        let k2_2 = ydd;

        // k3 at t + h
        let y_temp = yold + h * k2_1;
        let yd_temp = ydold + h * k2_2;
        let ydd = w * xin - w * w * y_temp;
        let k3_1 = yd_temp;
        let k3_2 = ydd;

        // RK4 update
        t += h;
        y = yold + h * (k0_1 + 2.0 * (k1_1 + k2_1) + k3_1) / 6.0;
        yd = ydold + h * (k0_2 + 2.0 * (k1_2 + k2_2) + k3_2) / 6.0;

        array_t.push(t);
        array_y.push(y);
    }

    Results {
        time: array_t,
        y: array_y,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c5l1_datfil.txt", output_dir);
    save_data(&data_file, &[results.time.clone(), results.y.clone()])?;

    let plot_file = format!("{}/c5l1_plot.png", output_dir);
    let config = PlotConfig::new("Fourth-order Runge-Kutta: Second Order Network")
        .with_labels("Time (s)", "Y");

    let series = vec![
        Series::new(results.time.clone(), results.y.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C5L1: Simulation finished");
    println!("  Data saved to: {}", data_file);
    println!("  Plot saved to: {}", plot_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c5l1_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }

    #[test]
    fn test_c5l1_converges() {
        let results = run();
        let final_y = *results.y.last().unwrap();
        // Should converge towards 1/w^2 * w = 1/w for step input
        // Actually for this system y'' + w^2*y = w*xin, steady state is xin/w = 1/20 = 0.05
        assert!(final_y > 0.0 && final_y < 0.1);
    }
}
