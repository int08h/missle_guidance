//! Chapter 26, Lesson 10: Frequency Response - Alternative Formulation
//!
//! Computes gain and phase frequency response for an alternative transfer function.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub freq: Vec<f64>,    // Frequency (rad/sec)
    pub gain: Vec<f64>,    // Gain (dB)
    pub phase: Vec<f64>,   // Phase (deg)
}

/// Run the C26L10 simulation
pub fn run() -> Results {
    let mut array_w = Vec::new();
    let mut array_gain = Vec::new();
    let mut array_phase = Vec::new();

    for i in 2..=160 {
        let w = 10.0_f64.powf(0.025 * (i as f64) - 1.0);
        let c1: f64 = 1.0;
        let c2: f64 = 0.0363;
        let top = 9000.0 * (c1.powi(2) + (c2 * w).powi(2)).sqrt();
        let bot = 2.0 * w * (1.0 + (w / 2.0).powi(2)).sqrt();
        let xmag = top / bot;
        let gain = 20.0 * xmag.log10();
        let phase = 57.3 * (w).atan2(29.1) - 90.0 - 57.3 * (w).atan2(2.0);

        array_w.push(w);
        array_gain.push(gain);
        array_phase.push(phase);
    }

    Results {
        freq: array_w,
        gain: array_gain,
        phase: array_phase,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c26l10_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.freq.clone(),
        results.gain.clone(),
        results.phase.clone(),
    ])?;

    let plot_file = format!("{}/c26l10_gain.png", output_dir);
    let config = PlotConfig::new("Frequency Response (Alternative) - Gain")
        .with_labels("Frequency (Rad/Sec)", "Gain (dB)");

    let series = vec![
        Series::new(results.freq.clone(), results.gain.clone())
            .with_color(plotters::prelude::BLUE)
            .with_label("Gain"),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C26L10: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c26l10_runs() {
        let results = run();
        assert!(!results.freq.is_empty());
    }
}
