//! Chapter 16, Lesson 1: Two-Stage Rocket Performance
//!
//! Calculates staging parameters and simulates velocity buildup
//! for a two-stage rocket with given delta-V requirements.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub vk: Vec<f64>,     // Velocity (kft/s)
    pub ag: Vec<f64>,     // Acceleration (G)
}

/// Run the C16L1 simulation
pub fn run() -> Results {
    let xisp1: f64 = 250.0;    // First stage specific impulse
    let xisp2: f64 = 250.0;    // Second stage specific impulse
    let xmf1: f64 = 0.85;      // First stage mass fraction
    let xmf2: f64 = 0.85;      // Second stage mass fraction
    let wpay: f64 = 100.0;     // Payload weight (lb)
    let delv: f64 = 20000.0;   // Total delta-V (ft/s)
    let delv1 = 0.3333 * delv; // First stage delta-V
    let delv2 = 0.6667 * delv; // Second stage delta-V
    let amax1: f64 = 10.0;     // Max acceleration stage 1 (G)
    let amax2: f64 = 10.0;     // Max acceleration stage 2 (G)
    let h: f64 = 0.01;

    // Stage 2 sizing
    let top2 = wpay * ((delv2 / (xisp2 * 32.2)).exp() - 1.0);
    let bot2 = 1.0 / xmf2 - ((1.0 - xmf2) / xmf2) * (delv2 / (xisp2 * 32.2)).exp();
    let wp2 = top2 / bot2;
    let ws2 = wp2 * (1.0 - xmf2) / xmf2;
    let wtot2 = wp2 + ws2 + wpay;
    let trst2 = amax2 * (wpay + ws2);
    let tb2 = xisp2 * wp2 / trst2;

    // Stage 1 sizing
    let top1 = wtot2 * ((delv1 / (xisp1 * 32.2)).exp() - 1.0);
    let bot1 = 1.0 / xmf1 - ((1.0 - xmf1) / xmf1) * (delv1 / (xisp1 * 32.2)).exp();
    let wp1 = top1 / bot1;
    let ws1 = wp1 * (1.0 - xmf1) / xmf1;
    let wtot = wp1 + ws1 + wtot2;
    let trst1 = amax1 * (wtot2 + ws1);
    let tb1 = xisp1 * wp1 / trst1;

    let mut t: f64 = 0.0;
    let mut s: f64 = 0.0;
    let mut v: f64 = 0.0;

    let mut array_t = Vec::new();
    let mut array_vk = Vec::new();
    let mut array_ag = Vec::new();

    while t <= tb1 + tb2 {
        let vold = v;

        // Compute weight and thrust
        let (wgt, trst) = if t < tb1 {
            let wgt = -wp1 * t / tb1 + wtot;
            (wgt, trst1)
        } else if t < tb1 + tb2 {
            let wgt = -wp2 * t / tb2 + wtot2 + wp2 * tb1 / tb2;
            (wgt, trst2)
        } else {
            (wpay, 0.0)
        };

        let a = 32.2 * trst / wgt;

        // Euler step
        v += h * a;
        t += h;

        // Second derivative for RK2
        let (wgt, trst) = if t < tb1 {
            let wgt = -wp1 * t / tb1 + wtot;
            (wgt, trst1)
        } else if t < tb1 + tb2 {
            let wgt = -wp2 * t / tb2 + wtot2 + wp2 * tb1 / tb2;
            (wgt, trst2)
        } else {
            (wpay, 0.0)
        };

        let a = 32.2 * trst / wgt;

        // RK2 averaging
        v = (vold + v) / 2.0 + 0.5 * h * a;

        s += h;

        if s >= 0.99999 {
            s = 0.0;
            let ag = a / 32.2;
            let vk = v / 1000.0;

            array_t.push(t);
            array_vk.push(vk);
            array_ag.push(ag);
        }
    }

    Results {
        time: array_t,
        vk: array_vk,
        ag: array_ag,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c16l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.vk.clone(),
        results.ag.clone(),
    ])?;

    let plot_file1 = format!("{}/c16l1_velocity.png", output_dir);
    let config = PlotConfig::new("Two-Stage Rocket - Velocity")
        .with_labels("Time (Sec)", "Velocity (Kft/Sec)");

    let series = vec![
        Series::new(results.time.clone(), results.vk.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file1, &config, &series).ok();

    let plot_file2 = format!("{}/c16l1_accel.png", output_dir);
    let config = PlotConfig::new("Two-Stage Rocket - Acceleration")
        .with_labels("Time (Sec)", "Acceleration (G)");

    let series = vec![
        Series::new(results.time.clone(), results.ag.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file2, &config, &series).ok();

    println!("C16L1: Simulation finished");
    println!("  Data saved to: {}", data_file);
    println!("  Final velocity: {:.2} kft/s", results.vk.last().unwrap_or(&0.0));

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c16l1_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }

    #[test]
    fn test_c16l1_achieves_deltav() {
        let results = run();
        let final_v = results.vk.last().unwrap_or(&0.0) * 1000.0;
        // Should achieve approximately 20 kft/s delta-V
        assert!(final_v > 19000.0 && final_v < 21000.0);
    }
}
