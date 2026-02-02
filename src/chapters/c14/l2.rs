//! Chapter 14, Lesson 2: Command Guidance
//!
//! 2D engagement simulation with command guidance (THETT-THETM error).

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub rt1k: Vec<f64>,
    pub rt2k: Vec<f64>,
    pub rm1k: Vec<f64>,
    pub rm2k: Vec<f64>,
    pub xncg: Vec<f64>,
}

/// Run the C14L2 simulation
pub fn run() -> Results {
    let vm: f64 = 3000.0;
    let vt: f64 = 1000.0;
    let _xnt: f64 = 0.0;
    let rm1ic: f64 = 0.0;
    let rm2ic: f64 = 1.0;
    let rt1ic: f64 = 40000.0;
    let rt2ic: f64 = 10000.0;
    let hedeg: f64 = 0.0;
    let xnp: f64 = 10.0;
    let ts: f64 = 0.1;

    let mut rt1 = rt1ic;
    let mut rt2 = rt2ic;
    let mut rm1 = rm1ic;
    let mut rm2 = rm2ic;
    let mut betat: f64 = 0.0;
    let mut vt1 = -vt * betat.cos();
    let mut vt2 = vt * betat.sin();
    let he = hedeg / 57.3;
    let mut xnc: f64;
    let mut t: f64 = 0.0;
    let mut s: f64 = 0.0;
    let mut h: f64;

    let mut rtm1 = rt1 - rm1;
    let mut rtm2 = rt2 - rm2;
    let mut rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
    let thett = rt2.atan2(rt1);
    let mut vm1 = vm * (thett + he).cos();
    let mut vm2 = vm * (thett + he).sin();
    let vtm1 = vt1 - vm1;
    let vtm2 = vt2 - vm2;
    let mut vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;

    let mut array_t = Vec::new();
    let mut array_rt1k = Vec::new();
    let mut array_rt2k = Vec::new();
    let mut array_rm1k = Vec::new();
    let mut array_rm2k = Vec::new();
    let mut array_xncg = Vec::new();

    let mut am1: f64;
    let mut am2: f64;
    let xnt: f64 = 0.0;

    while vc >= 0.0 {
        if rtm < 1000.0 {
            h = 0.0002;
        } else {
            h = 0.01;
        }

        let betat_old = betat;
        let rt1_old = rt1;
        let rt2_old = rt2;
        let rm1_old = rm1;
        let rm2_old = rm2;
        let vm1_old = vm1;
        let vm2_old = vm2;

        // First derivative evaluation
        let thett = rt2.atan2(rt1);
        let thetm = rm2.atan2(rm1);
        let _rt = (rt1 * rt1 + rt2 * rt2).sqrt();
        let rm = (rm1 * rm1 + rm2 * rm2).sqrt();
        rtm1 = rt1 - rm1;
        rtm2 = rt2 - rm2;
        rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
        let vtm1 = vt1 - vm1;
        let vtm2 = vt2 - vm2;
        let xlam = rtm2.atan2(rtm1);
        let _xlamd = (rtm1 * vtm2 - rtm2 * vtm1) / (rtm * rtm);
        xnc = xnp * rm * (thett - thetm);
        am1 = -xnc * xlam.sin();
        am2 = xnc * xlam.cos();
        vt1 = -vt * betat.cos();
        vt2 = vt * betat.sin();
        let betatd = xnt / vt;

        // Euler step
        betat += h * betatd;
        rt1 += h * vt1;
        rt2 += h * vt2;
        rm1 += h * vm1;
        rm2 += h * vm2;
        vm1 += h * am1;
        vm2 += h * am2;
        t += h;

        // Second derivative
        let thett = rt2.atan2(rt1);
        let thetm = rm2.atan2(rm1);
        let _rt = (rt1 * rt1 + rt2 * rt2).sqrt();
        let rm = (rm1 * rm1 + rm2 * rm2).sqrt();
        rtm1 = rt1 - rm1;
        rtm2 = rt2 - rm2;
        rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
        let vtm1 = vt1 - vm1;
        let vtm2 = vt2 - vm2;
        vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;
        let xlam = rtm2.atan2(rtm1);
        let _xlamd = (rtm1 * vtm2 - rtm2 * vtm1) / (rtm * rtm);
        xnc = xnp * rm * (thett - thetm);
        am1 = -xnc * xlam.sin();
        am2 = xnc * xlam.cos();
        vt1 = -vt * betat.cos();
        vt2 = vt * betat.sin();
        let betatd = xnt / vt;

        // RK2 averaging
        betat = 0.5 * (betat_old + betat + h * betatd);
        rt1 = 0.5 * (rt1_old + rt1 + h * vt1);
        rt2 = 0.5 * (rt2_old + rt2 + h * vt2);
        rm1 = 0.5 * (rm1_old + rm1 + h * vm1);
        rm2 = 0.5 * (rm2_old + rm2 + h * vm2);
        vm1 = 0.5 * (vm1_old + vm1 + h * am1);
        vm2 = 0.5 * (vm2_old + vm2 + h * am2);

        s += h;
        if s >= (ts - 1e-5) {
            s = 0.0;
            let rt1k = rt1 / 1000.0;
            let rt2k = rt2 / 1000.0;
            let rm1k = rm1 / 1000.0;
            let rm2k = rm2 / 1000.0;

            array_t.push(t);
            array_rt1k.push(rt1k);
            array_rt2k.push(rt2k);
            array_rm1k.push(rm1k);
            array_rm2k.push(rm2k);
            array_xncg.push(xnc / 32.2);
        }
    }

    println!("RTM = {:.6}", rtm);

    Results {
        time: array_t,
        rt1k: array_rt1k,
        rt2k: array_rt2k,
        rm1k: array_rm1k,
        rm2k: array_rm2k,
        xncg: array_xncg,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c14l2_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.rt1k.clone(),
        results.rt2k.clone(),
        results.rm1k.clone(),
        results.rm2k.clone(),
        results.xncg.clone(),
    ])?;

    // Trajectory plot
    let plot_file = format!("{}/c14l2_trajectory.png", output_dir);
    let config = PlotConfig::new("Trajectory")
        .with_labels("Downrange (Kft)", "Altitude (Kft)");
    let series = vec![
        Series::new(results.rt1k.clone(), results.rt2k.clone())
            .with_label("Target")
            .with_color(plotters::prelude::BLUE),
        Series::new(results.rm1k.clone(), results.rm2k.clone())
            .with_label("Missile")
            .with_color(plotters::prelude::RED),
    ];
    line_plot(&plot_file, &config, &series).ok();

    // Acceleration plot
    let plot_file2 = format!("{}/c14l2_acceleration.png", output_dir);
    let config2 = PlotConfig::new("Missile Acceleration")
        .with_labels("Time (S)", "Acceleration (G)");
    let series2 = vec![
        Series::new(results.time.clone(), results.xncg.clone())
            .with_color(plotters::prelude::BLUE),
    ];
    line_plot(&plot_file2, &config2, &series2).ok();

    println!("C14L2: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c14l2_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
