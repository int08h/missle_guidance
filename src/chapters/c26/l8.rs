//! Chapter 26, Lesson 8: Frequency Response with Compensator
//!
//! Computes gain frequency response with actuator compensator.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub freq: Vec<f64>,    // Frequency (rad/sec)
    pub gain: Vec<f64>,    // Gain (dB)
}

/// Run the C26L8 simulation
pub fn run() -> Results {
    let mut array_w = Vec::new();
    let mut array_gain = Vec::new();

    for i in 2..=160 {
        let w = 10.0_f64.powf(0.025 * (i as f64) - 1.0);
        let top1 = 13500.0 * (1.0 + (w / 29.1).powi(2)).sqrt();
        let top2 = ((1.0 - w * w / (173.0 * 173.0)).powi(2) + (2.0 * 0.74 * w / 173.0).powi(2)).sqrt();
        let bot1 = w * (1.0 + (w / 2.0).powi(2)).sqrt();
        let bot2 = ((1.0 - w * w / (100.0 * 100.0)).powi(2) + (2.0 * 0.65 * w / 100.0).powi(2)).sqrt();
        let xmag = top1 * top2 / (bot1 * bot2);
        let gain = 20.0 * xmag.log10();

        array_w.push(w);
        array_gain.push(gain);
    }

    Results {
        freq: array_w,
        gain: array_gain,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c26l8_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.freq.clone(),
        results.gain.clone(),
    ])?;

    let plot_file = format!("{}/c26l8_gain.png", output_dir);
    let config = PlotConfig::new("Frequency Response with Compensator - Gain")
        .with_labels("Frequency (Rad/Sec)", "Gain (dB)");

    let series = vec![
        Series::new(results.freq.clone(), results.gain.clone())
            .with_color(plotters::prelude::BLUE)
            .with_label("Gain"),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C26L8: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c26l8_runs() {
        let results = run();
        assert!(!results.freq.is_empty());
    }
}
