//! Chapter 11, Lesson 1: Reentry Target Trajectory
//!
//! Simulates a ballistic reentry vehicle trajectory with atmospheric drag.
//! Initial conditions: 100 kft altitude, 6000 ft/s velocity, 45 deg flight path.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub rt1k: Vec<f64>,    // Downrange (Kft)
    pub rt2k: Vec<f64>,    // Altitude (Kft)
    pub atg: Vec<f64>,     // Acceleration (G)
    pub vt: Vec<f64>,      // Velocity (ft/s)
}

/// Run the C11L1 simulation
pub fn run() -> Results {
    let h: f64 = 0.01;
    let beta: f64 = 500.0;     // Ballistic coefficient
    let vt_init: f64 = 6000.0;
    let gamtdeg: f64 = 45.0;   // Flight path angle (deg)

    let mut t: f64 = 0.0;
    let mut s: f64 = 0.0;

    // Initial conditions
    let mut rt1: f64 = 0.0;
    let mut rt2: f64 = 100000.0;
    let mut vt1 = vt_init * (gamtdeg / 57.3).cos();
    let mut vt2 = -vt_init * (gamtdeg / 57.3).sin();

    let mut array_t = Vec::new();
    let mut array_rt1k = Vec::new();
    let mut array_rt2k = Vec::new();
    let mut array_atg = Vec::new();
    let mut array_vt = Vec::new();

    while rt2 >= 0.0 {
        let rt1old = rt1;
        let rt2old = rt2;
        let vt1old = vt1;
        let vt2old = vt2;

        // First derivative evaluation
        let rho = if rt2 <= 30000.0 {
            0.002378 * (-rt2 / 30000.0).exp()
        } else {
            0.0034 * (-rt2 / 22000.0).exp()
        };

        let vt = (vt1 * vt1 + vt2 * vt2).sqrt();
        let q = 0.5 * rho * vt * vt;
        let gamt = (-vt2).atan2(vt1);
        let at1 = -32.2 * q * gamt.cos() / beta;
        let at2 = -32.2 + 32.2 * q * gamt.sin() / beta;

        // Euler step
        rt1 += h * vt1;
        rt2 += h * vt2;
        vt1 += h * at1;
        vt2 += h * at2;
        t += h;

        // Second derivative for RK2
        let rho = if rt2 <= 30000.0 {
            0.002378 * (-rt2 / 30000.0).exp()
        } else {
            0.0034 * (-rt2 / 22000.0).exp()
        };

        let vt = (vt1 * vt1 + vt2 * vt2).sqrt();
        let q = 0.5 * rho * vt * vt;
        let gamt = (-vt2).atan2(vt1);
        let at1 = -32.2 * q * gamt.cos() / beta;
        let at2 = -32.2 + 32.2 * q * gamt.sin() / beta;

        // RK2 averaging
        rt1 = 0.5 * (rt1old + rt1 + h * vt1);
        rt2 = 0.5 * (rt2old + rt2 + h * vt2);
        vt1 = 0.5 * (vt1old + vt1 + h * at1);
        vt2 = 0.5 * (vt2old + vt2 + h * at2);

        s += h;

        if s >= 0.09999 {
            s = 0.0;
            let atg = (at1 * at1 + at2 * at2).sqrt() / 32.2;
            let rt1k = rt1 / 1000.0;
            let rt2k = rt2 / 1000.0;
            let vt = (vt1 * vt1 + vt2 * vt2).sqrt();

            array_t.push(t);
            array_rt1k.push(rt1k);
            array_rt2k.push(rt2k);
            array_atg.push(atg);
            array_vt.push(vt);
        }
    }

    Results {
        time: array_t,
        rt1k: array_rt1k,
        rt2k: array_rt2k,
        atg: array_atg,
        vt: array_vt,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c11l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.rt1k.clone(),
        results.rt2k.clone(),
        results.atg.clone(),
        results.vt.clone(),
    ])?;

    let plot_file = format!("{}/c11l1_trajectory.png", output_dir);
    let config = PlotConfig::new("Reentry Target Trajectory")
        .with_labels("Downrange (Kft)", "Altitude (Kft)");

    let series = vec![
        Series::new(results.rt1k.clone(), results.rt2k.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C11L1: Simulation finished");
    println!("  Data saved to: {}", data_file);
    println!("  Plot saved to: {}", plot_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c11l1_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }

    #[test]
    fn test_c11l1_trajectory_ends_at_ground() {
        let results = run();
        let final_alt = *results.rt2k.last().unwrap();
        assert!(final_alt <= 0.1);  // Within 100 ft of ground
    }
}
