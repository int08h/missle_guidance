//! Chapter 10, Lesson 1: Ballistic Trajectory with Drag
//!
//! 2D point mass ballistic trajectory simulation with altitude-dependent
//! atmospheric density and aerodynamic drag.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub rm1k: Vec<f64>,  // Downrange (Kft)
    pub rm2k: Vec<f64>,  // Altitude (Kft)
}

/// Run the C10L1 simulation
pub fn run() -> Results {
    let h: f64 = 0.01;
    let vm: f64 = 3000.0;
    let beta: f64 = 1000.0;  // Ballistic coefficient
    let gamdeg: f64 = 45.0;  // Launch angle

    let mut t: f64 = 0.0;
    let mut s: f64 = 0.0;

    // Initial conditions
    let mut vm1 = vm * (gamdeg / 57.3).cos();
    let mut vm2 = vm * (gamdeg / 57.3).sin();
    let mut rm1: f64 = 0.0;
    let mut rm2: f64 = 0.0;

    let mut array_t = Vec::new();
    let mut array_rm1k = Vec::new();
    let mut array_rm2k = Vec::new();

    // Integrate until impact (rm2 <= 0 after launch)
    while !(t > 0.0 && rm2 <= 0.0) {
        let rm1old = rm1;
        let rm2old = rm2;
        let vm1old = vm1;
        let vm2old = vm2;

        // First derivative evaluation
        let rho = if rm2 < 30000.0 {
            0.002378 * (-rm2 / 30000.0).exp()
        } else {
            0.0034 * (-rm2 / 22000.0).exp()
        };

        let vm_mag = (vm1 * vm1 + vm2 * vm2).sqrt();
        let q = 0.5 * rho * vm_mag * vm_mag;
        let gam = vm2.atan2(vm1);
        let drag = q * 32.2 / beta;
        let am1 = -drag * gam.cos();
        let am2 = -32.2 - drag * gam.sin();

        // Euler step
        rm1 += h * vm1;
        rm2 += h * vm2;
        vm1 += h * am1;
        vm2 += h * am2;
        t += h;

        // Second derivative for RK2
        let rho = if rm2 < 30000.0 {
            0.002378 * (-rm2 / 30000.0).exp()
        } else {
            0.0034 * (-rm2 / 22000.0).exp()
        };

        let vm_mag = (vm1 * vm1 + vm2 * vm2).sqrt();
        let q = 0.5 * rho * vm_mag * vm_mag;
        let gam = vm2.atan2(vm1);
        let drag = q * 32.2 / beta;
        let am1 = -drag * gam.cos();
        let am2 = -32.2 - drag * gam.sin();

        // RK2 averaging
        rm1 = 0.5 * (rm1old + rm1 + h * vm1);
        rm2 = 0.5 * (rm2old + rm2 + h * vm2);
        vm1 = 0.5 * (vm1old + vm1 + h * am1);
        vm2 = 0.5 * (vm2old + vm2 + h * am2);

        s += h;

        if s >= 0.099999 {
            s = 0.0;
            let rm1k = rm1 / 1000.0;
            let rm2k = rm2 / 1000.0;
            array_t.push(t);
            array_rm1k.push(rm1k);
            array_rm2k.push(rm2k);
        }
    }

    Results {
        time: array_t,
        rm1k: array_rm1k,
        rm2k: array_rm2k,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c10l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.rm1k.clone(),
        results.rm2k.clone(),
    ])?;

    let plot_file = format!("{}/c10l1_plot.png", output_dir);
    let config = PlotConfig::new("Ballistic Trajectory with Drag")
        .with_labels("Downrange (Kft)", "Altitude (Kft)");

    let series = vec![
        Series::new(results.rm1k.clone(), results.rm2k.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C10L1: Simulation finished");
    println!("  Data saved to: {}", data_file);
    println!("  Plot saved to: {}", plot_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c10l1_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }

    #[test]
    fn test_c10l1_impacts() {
        let results = run();
        // Trajectory should end near ground level
        let final_alt = *results.rm2k.last().unwrap();
        assert!(final_alt <= 0.1);  // Within 100 ft of ground
    }

    #[test]
    fn test_c10l1_positive_range() {
        let results = run();
        // All downrange values should be positive
        for rm1k in &results.rm1k {
            assert!(*rm1k >= 0.0);
        }
    }
}
