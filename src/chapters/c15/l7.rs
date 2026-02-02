//! Chapter 15, Lesson 7: Velocity and Flight Time vs Flight Path Angle
//!
//! Parametric study of velocity and flight time for 10000 km range.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub gamdeg: Vec<f64>,
    pub vkm: Vec<f64>,
    pub tf: Vec<f64>,
    pub xlam: Vec<f64>,
}

/// Run the C15L7 simulation
pub fn run() -> Results {
    let gm: f64 = 1.4077e16;
    let a: f64 = 2.0926e7;
    let distkm: f64 = 10000.0;
    let phi = distkm * 3280.0 / a;

    let mut array_gamdeg = Vec::new();
    let mut array_vkm = Vec::new();
    let mut array_tf = Vec::new();
    let mut array_xlam = Vec::new();

    // Sweep from 160 to 180 degrees
    for gamdeg_val in 160..=180 {
        let gamdeg = gamdeg_val as f64;
        let gam = gamdeg / 57.3;

        let top = gm * (1.0 - phi.cos());
        let temp = a * gam.cos() / a - (phi + gam).cos();
        let bot = a * gam.cos() * temp;
        let v = (top / bot).sqrt();
        let vkm = v / 3280.0;
        let xlam = a * v * v / gm;

        let top1 = gam.tan() * (1.0 - phi.cos()) + (1.0 - xlam) * phi.sin();
        let bot1p = (1.0 - phi.cos()) / (xlam * gam.cos() * gam.cos());
        let bot1 = (2.0 - xlam) * (bot1p + (gam + phi).cos() / gam.cos());
        let top2 = 2.0 * gam.cos();
        let bot2 = xlam * (2.0 / xlam - 1.0).powf(1.5);
        let top3 = (2.0 / xlam - 1.0).sqrt();
        let bot3 = gam.cos() / (phi / 2.0).tan() - gam.sin();
        let temp = (top2 / bot2) * top3.atan2(bot3);
        let tf = a * (top1 / bot1 + temp) / (v * gam.cos());

        array_gamdeg.push(gamdeg);
        array_vkm.push(vkm);
        array_tf.push(tf);
        array_xlam.push(xlam);
    }

    Results {
        gamdeg: array_gamdeg,
        vkm: array_vkm,
        tf: array_tf,
        xlam: array_xlam,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c15l7_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.gamdeg.clone(),
        results.vkm.clone(),
        results.tf.clone(),
        results.xlam.clone(),
    ])?;

    let plot_file = format!("{}/c15l7_velocity.png", output_dir);
    let config = PlotConfig::new("Velocity vs Flight Path Angle")
        .with_labels("Flight Path Angle (deg)", "Velocity (km/s)");
    let series = vec![
        Series::new(results.gamdeg.clone(), results.vkm.clone())
            .with_color(plotters::prelude::BLUE),
    ];
    line_plot(&plot_file, &config, &series).ok();

    let plot_file2 = format!("{}/c15l7_flight_time.png", output_dir);
    let config2 = PlotConfig::new("Flight Time vs Flight Path Angle")
        .with_labels("Flight Path Angle (deg)", "Flight Time (s)");
    let series2 = vec![
        Series::new(results.gamdeg.clone(), results.tf.clone())
            .with_color(plotters::prelude::BLUE),
    ];
    line_plot(&plot_file2, &config2, &series2).ok();

    println!("C15L7: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c15l7_runs() {
        let results = run();
        assert!(!results.gamdeg.is_empty());
    }
}
