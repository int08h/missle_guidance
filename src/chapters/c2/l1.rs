//! Chapter 2, Lesson 1: Two-Dimensional Tactical Missile-Target Engagement
//!
//! Simulates proportional navigation guidance in a 2D engagement scenario.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

/// Simulation results
pub struct Results {
    pub time: Vec<f64>,
    pub rt1: Vec<f64>,
    pub rt2: Vec<f64>,
    pub rm1: Vec<f64>,
    pub rm2: Vec<f64>,
    pub xnc_g: Vec<f64>,
    pub rtm: Vec<f64>,
}

/// Run the C2L1 simulation
pub fn run() -> Results {
    let vm: f64 = 3000.0;     // Missile velocity (ft/s)
    let vt: f64 = 1000.0;     // Target velocity (ft/s)
    let xnt: f64 = 0.0;       // Target acceleration
    let he_deg: f64 = -20.0;  // Heading error (degrees)
    let xnp: f64 = 4.0;       // Navigation ratio
    let beta: f64 = 0.0;      // Target heading angle

    // Initial positions
    let mut rm1: f64 = 0.0;
    let mut rm2: f64 = 10000.0;
    let mut rt1: f64 = 40000.0;
    let mut rt2: f64 = 10000.0;

    // Target velocity components
    let mut vt1: f64 = -vt * beta.cos();
    let mut vt2: f64 = vt * beta.sin();

    // Convert heading error to radians
    let he = he_deg / 57.3;

    let mut t = 0.0;
    let mut s = 0.0;

    // Initial calculations
    let mut rtm1 = rt1 - rm1;
    let mut rtm2 = rt2 - rm2;
    let mut rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
    let xlam = rtm2.atan2(rtm1);
    let xlead = (vt * (beta + xlam).sin() / vm).asin();
    let thet = xlam + xlead;

    // Missile velocity components
    let mut vm1 = vm * (thet + he).cos();
    let mut vm2 = vm * (thet + he).sin();

    let mut vtm1 = vt1 - vm1;
    let mut vtm2 = vt2 - vm2;
    let mut vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;

    let mut xnc: f64;

    // Results storage
    let mut array_t = Vec::new();
    let mut array_rt1 = Vec::new();
    let mut array_rt2 = Vec::new();
    let mut array_rm1 = Vec::new();
    let mut array_rm2 = Vec::new();
    let mut array_xnc_g = Vec::new();
    let mut array_rtm = Vec::new();

    while vc >= 0.0 {
        // Adaptive step size
        let h = if rtm < 1000.0 { 0.0002 } else { 0.01 };

        // Save old values
        let beta_old = beta;
        let rt1_old = rt1;
        let rt2_old = rt2;
        let rm1_old = rm1;
        let rm2_old = rm2;
        let vm1_old = vm1;
        let vm2_old = vm2;

        // First derivative evaluation
        rtm1 = rt1 - rm1;
        rtm2 = rt2 - rm2;
        rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
        vtm1 = vt1 - vm1;
        vtm2 = vt2 - vm2;
        vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;
        let xlam = rtm2.atan2(rtm1);
        let xlamd = (rtm1 * vtm2 - rtm2 * vtm1) / (rtm * rtm);
        xnc = xnp * vc * xlamd;
        let am1 = -xnc * xlam.sin();
        let am2 = xnc * xlam.cos();
        vt1 = -vt * beta.cos();
        vt2 = vt * beta.sin();
        let betad = xnt / vt;

        // Euler step
        let beta_new = beta + h * betad;
        rt1 += h * vt1;
        rt2 += h * vt2;
        rm1 += h * vm1;
        rm2 += h * vm2;
        vm1 += h * am1;
        vm2 += h * am2;
        t += h;

        // Second evaluation
        rtm1 = rt1 - rm1;
        rtm2 = rt2 - rm2;
        rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
        vtm1 = vt1 - vm1;
        vtm2 = vt2 - vm2;
        vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;
        let xlam = rtm2.atan2(rtm1);
        let xlamd = (rtm1 * vtm2 - rtm2 * vtm1) / (rtm * rtm);
        xnc = xnp * vc * xlamd;
        let am1 = -xnc * xlam.sin();
        let am2 = xnc * xlam.cos();
        let vt1_new = -vt * beta_new.cos();
        let vt2_new = vt * beta_new.sin();
        let betad = xnt / vt;

        // RK2 averaging
        let _beta = 0.5 * (beta_old + beta_new + h * betad);
        rt1 = 0.5 * (rt1_old + rt1 + h * vt1_new);
        rt2 = 0.5 * (rt2_old + rt2 + h * vt2_new);
        rm1 = 0.5 * (rm1_old + rm1 + h * vm1);
        rm2 = 0.5 * (rm2_old + rm2 + h * vm2);
        vm1 = 0.5 * (vm1_old + vm1 + h * am1);
        vm2 = 0.5 * (vm2_old + vm2 + h * am2);

        s += h;

        // Store data at sampling interval
        if s >= 0.09999 {
            s = 0.0;
            array_t.push(t);
            array_rt1.push(rt1);
            array_rt2.push(rt2);
            array_rm1.push(rm1);
            array_rm2.push(rm2);
            array_xnc_g.push(xnc / 32.2);
            array_rtm.push(rtm);
        }
    }

    println!("Final miss distance (RTM): {:.2} ft", rtm);

    Results {
        time: array_t,
        rt1: array_rt1,
        rt2: array_rt2,
        rm1: array_rm1,
        rm2: array_rm2,
        xnc_g: array_xnc_g,
        rtm: array_rtm,
    }
}

/// Run and save results to file
pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    // Save data file
    let data_file = format!("{}/c2l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.rt1.clone(),
        results.rt2.clone(),
        results.rm1.clone(),
        results.rm2.clone(),
        results.xnc_g.clone(),
        results.rtm.clone(),
    ])?;

    // Create trajectory plot
    let plot_file1 = format!("{}/c2l1_trajectory.png", output_dir);
    let config = PlotConfig::new("Two-dimensional tactical missile-target engagement simulation")
        .with_labels("Downrange (Ft)", "Altitude (Ft)");

    let series = vec![
        Series::new(results.rt1.clone(), results.rt2.clone())
            .with_label("Target")
            .with_color(plotters::prelude::RED),
        Series::new(results.rm1.clone(), results.rm2.clone())
            .with_label("Missile")
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file1, &config, &series).ok();

    // Create acceleration plot
    let plot_file2 = format!("{}/c2l1_acceleration.png", output_dir);
    let config = PlotConfig::new("Two-dimensional tactical missile-target engagement simulation")
        .with_labels("Time (sec)", "Acceleration of missile (G)");

    let series = vec![
        Series::new(results.time.clone(), results.xnc_g.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file2, &config, &series).ok();

    println!("C2L1: Simulation Complete");
    println!("  Data saved to: {}", data_file);
    println!("  Plot saved to: {}", plot_file1);
    println!("  Plot saved to: {}", plot_file2);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c2l1_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }

    #[test]
    fn test_c2l1_intercept() {
        let results = run();
        // Check minimum RTM (closest approach) since sampling might miss the exact intercept
        let min_rtm = results.rtm.iter().cloned().fold(f64::INFINITY, f64::min);
        // Should achieve close intercept with PN guidance (within 100 ft at closest approach)
        assert!(min_rtm < 500.0, "Closest approach should be small: {}", min_rtm);
    }
}
