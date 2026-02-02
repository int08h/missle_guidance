//! Chapter 22, Lesson 3: Frequency Response Analysis (Bode Plot)
//!
//! Computes frequency response (gain and phase) for the rate gyro
//! flight control system using analytical transfer function.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub freq: Vec<f64>,    // Frequency (rad/sec)
    pub gain: Vec<f64>,    // Gain (dB)
    pub phase: Vec<f64>,   // Phase (deg)
}

/// Run the C22L3 simulation
pub fn run() -> Results {
    let zact: f64 = 0.7;
    let wact: f64 = 150.0;
    let k3: f64 = -1.89;
    let ta: f64 = 0.457;
    let zaf: f64 = 0.058;
    let waf: f64 = 25.3;
    let kr: f64 = 0.1;

    let mut array_w = Vec::new();
    let mut array_gain = Vec::new();
    let mut array_phase = Vec::new();

    for i in 2..=160 {
        let w = 10.0_f64.powf(0.025 * (i as f64) - 1.0);

        let xmag1 = (1.0 + (w * ta).powi(2)).sqrt();
        let xmag2 = ((1.0 - (w / waf).powi(2)).powi(2) + (2.0 * zaf * w / waf).powi(2)).sqrt();
        let xmag3 = ((1.0 - (w / wact).powi(2)).powi(2) + (2.0 * zact * w / wact).powi(2)).sqrt();

        let gain = 20.0 * (-k3 * kr * xmag1 / (xmag2 * xmag3)).log10();

        let phase1 = 57.3 * (w * ta).atan2(1.0);
        let phase2 = 57.3 * (2.0 * zaf * w / waf).atan2(1.0 - (w / waf).powi(2));
        let phase3 = 57.3 * (2.0 * zact * w / wact).atan2(1.0 - (w / wact).powi(2));

        let phase = phase1 - phase2 - phase3;

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

    let data_file = format!("{}/c22l3_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.freq.clone(),
        results.gain.clone(),
        results.phase.clone(),
    ])?;

    let plot_file = format!("{}/c22l3_gain.png", output_dir);
    let config = PlotConfig::new("Frequency Response - Gain")
        .with_labels("Frequency (Rad/Sec)", "Gain (Db)");

    let series = vec![
        Series::new(results.freq.clone(), results.gain.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    let plot_file2 = format!("{}/c22l3_phase.png", output_dir);
    let config2 = PlotConfig::new("Frequency Response - Phase")
        .with_labels("Frequency (Rad/Sec)", "Phase (Deg)");

    let series2 = vec![
        Series::new(results.freq.clone(), results.phase.clone())
            .with_color(plotters::prelude::RED),
    ];

    line_plot(&plot_file2, &config2, &series2).ok();

    println!("C22L3: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c22l3_runs() {
        let results = run();
        assert!(!results.freq.is_empty());
        assert_eq!(results.freq.len(), 159);
    }
}
