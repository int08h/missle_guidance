//! Chapter 26, Lesson 9: Frequency Response with Gyro Compensation
//!
//! Computes gain frequency response with actuator and gyro compensators.
//!
//! MATLAB BUG: The output line includes `ArrayPHASE` but it is never computed
//! in the code. We output zeros for the phase column to match MATLAB's
//! uninitialized array behavior.
//! Output: 3 columns [W, GAIN, PHASE] (phase is always 0)

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub freq: Vec<f64>,    // Frequency (rad/sec)
    pub gain: Vec<f64>,    // Gain (dB)
    pub phase: Vec<f64>,   // Phase (deg) - not computed in MATLAB but in output
}

/// Run the C26L9 simulation
pub fn run() -> Results {
    let mut array_w = Vec::new();
    let mut array_gain = Vec::new();
    let mut array_phase = Vec::new();

    for i in 2..=160 {
        let w = 10.0_f64.powf(0.025 * (i as f64) - 1.0);
        let top1 = 13500.0 * (1.0 + (w / 29.0).powi(2)).sqrt();
        let top2 = ((1.0 - w * w / (160.0 * 160.0)).powi(2) + (2.0 * 0.83 * w / 160.0).powi(2)).sqrt();
        let top3 = ((1.0 - w * w / (216.0 * 216.0)).powi(2) + (2.0 * 0.45 * w / 216.0).powi(2)).sqrt();
        let bot1 = w * (1.0 + (w / 2.0).powi(2)).sqrt();
        let bot2 = ((1.0 - w * w / (100.0 * 100.0)).powi(2) + (2.0 * 0.65 * w / 100.0).powi(2)).sqrt();
        let bot3 = ((1.0 - w * w / (200.0 * 200.0)).powi(2) + (2.0 * 0.5 * w / 200.0).powi(2)).sqrt();
        let xmag = top1 * top2 * top3 / (bot1 * bot2 * bot3);
        let gain = 20.0 * xmag.log10();
        // MATLAB computes phase but ArrayPHASE is never assigned in C26L9
        // The output line includes ArrayPHASE but it would be empty/zeros
        // Looking at the code, phase is not computed - keep as 0
        let phase = 0.0;

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

    // Note: MATLAB outputs [ArrayW', ArrayGAIN', ArrayPHASE'] but ArrayPHASE is never set
    // We match MATLAB output format
    let data_file = format!("{}/c26l9_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.freq.clone(),
        results.gain.clone(),
        results.phase.clone(),
    ])?;

    let plot_file = format!("{}/c26l9_gain.png", output_dir);
    let config = PlotConfig::new("Frequency Response with Gyro Compensation - Gain")
        .with_labels("Frequency (Rad/Sec)", "Gain (dB)");

    let series = vec![
        Series::new(results.freq.clone(), results.gain.clone())
            .with_color(plotters::prelude::BLUE)
            .with_label("Gain"),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C26L9: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c26l9_runs() {
        let results = run();
        assert!(!results.freq.is_empty());
    }
}
