//! Chapter 15, Lesson 6: Flight Time vs Distance
//!
//! Analytical computation of flight time for minimum energy trajectory.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub distkm: Vec<f64>,
    pub tf: Vec<f64>,
}

/// Run the C15L6 simulation
pub fn run() -> Results {
    let gm: f64 = 1.4077e16;
    let a: f64 = 2.0926e7;

    let mut array_distkm = Vec::new();
    let mut array_tf = Vec::new();

    let mut distkm: f64 = 100.0;
    while distkm <= 20000.0 {
        let phi = distkm * 3280.0 / a;
        let gam = std::f64::consts::PI / 4.0 - phi / 4.0;

        let top = gm * (1.0 - phi.cos());
        let temp = a * gam.cos() / a - (phi + gam).cos();
        let bot = a * gam.cos() * temp;
        let v = (top / bot).sqrt();
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

        array_distkm.push(distkm);
        array_tf.push(tf);

        distkm += 100.0;
    }

    Results {
        distkm: array_distkm,
        tf: array_tf,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c15l6_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.distkm.clone(),
        results.tf.clone(),
    ])?;

    let plot_file = format!("{}/c15l6_flight_time.png", output_dir);
    let config = PlotConfig::new("Flight Time vs Distance")
        .with_labels("Distance (km)", "Flight Time (s)");
    let series = vec![
        Series::new(results.distkm.clone(), results.tf.clone())
            .with_color(plotters::prelude::BLUE),
    ];
    line_plot(&plot_file, &config, &series).ok();

    println!("C15L6: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c15l6_runs() {
        let results = run();
        assert!(!results.distkm.is_empty());
    }
}
