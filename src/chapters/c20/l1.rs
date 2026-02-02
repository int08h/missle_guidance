//! Chapter 20, Lesson 1: Target Jink/Displacement
//!
//! Simulates PN guidance response to a target that suddenly
//! displaces (jinks) perpendicular to the line of sight.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub rm1: Vec<f64>,
    pub rm2: Vec<f64>,
    pub rt1: Vec<f64>,
    pub rt2: Vec<f64>,
    pub xncg: Vec<f64>,
}

/// Run the C20L1 simulation
pub fn run() -> Results {
    let xnp: f64 = 3.0;
    let displace: f64 = 200.0;   // Target displacement (ft)
    let thom: f64 = 1.0;         // Time to home at jink
    let vm: f64 = 3000.0;
    let vt: f64 = 1000.0;

    // Initial conditions
    let mut rm1: f64 = 0.0;
    let mut rm2: f64 = 1000.0;
    let mut rt1: f64 = 20000.0;
    let mut rt2: f64 = 1000.0;

    let mut qswitch = false;
    let vt1 = -vt;
    let vt2: f64 = 0.0;

    let mut t: f64 = 0.0;
    let mut s: f64 = 0.0;

    let mut rtm1 = rt1 - rm1;
    let mut rtm2 = rt2 - rm2;
    let mut rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();

    let mut vm1 = vm;
    let mut vm2: f64 = 0.0;

    let vtm1 = vt1 - vm1;
    let vtm2 = vt2 - vm2;
    let mut vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;
    let mut tgo: f64;

    let mut h: f64 = 0.005;

    let mut array_t = Vec::new();
    let mut array_rm1 = Vec::new();
    let mut array_rm2 = Vec::new();
    let mut array_rt1 = Vec::new();
    let mut array_rt2 = Vec::new();
    let mut array_xncg = Vec::new();

    while vc > 0.0 {
        tgo = rtm / vc;
        if tgo < 0.3 {
            h = 0.00005;
        }

        // Target jink at specified time-to-go
        if tgo <= thom && !qswitch {
            qswitch = true;
            rt2 += displace;
        }

        let rt1old = rt1;
        let rt2old = rt2;
        let rm1old = rm1;
        let rm2old = rm2;
        let vm1old = vm1;
        let vm2old = vm2;

        // First derivative evaluation
        rtm1 = rt1 - rm1;
        rtm2 = rt2 - rm2;
        rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();

        let vtm1 = vt1 - vm1;
        let vtm2 = vt2 - vm2;
        vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;

        let xlam = rtm2.atan2(rtm1);
        let xlamd = (rtm1 * vtm2 - rtm2 * vtm1) / (rtm * rtm);
        let xnc = xnp * vc * xlamd;
        let am1 = -xnc * xlam.sin();
        let am2 = xnc * xlam.cos();

        // Euler step
        rt1 += h * vt1;
        rt2 += h * vt2;
        rm1 += h * vm1;
        rm2 += h * vm2;
        vm1 += h * am1;
        vm2 += h * am2;
        t += h;

        // Second derivative for RK2
        rtm1 = rt1 - rm1;
        rtm2 = rt2 - rm2;
        rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();

        let vtm1 = vt1 - vm1;
        let vtm2 = vt2 - vm2;
        vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;

        let xlam = rtm2.atan2(rtm1);
        let xlamd = (rtm1 * vtm2 - rtm2 * vtm1) / (rtm * rtm);
        let xnc = xnp * vc * xlamd;
        let am1 = -xnc * xlam.sin();
        let am2 = xnc * xlam.cos();

        // RK2 averaging
        rt1 = 0.5 * (rt1old + rt1 + h * vt1);
        rt2 = 0.5 * (rt2old + rt2 + h * vt2);
        rm1 = 0.5 * (rm1old + rm1 + h * vm1);
        rm2 = 0.5 * (rm2old + rm2 + h * vm2);
        vm1 = 0.5 * (vm1old + vm1 + h * am1);
        vm2 = 0.5 * (vm2old + vm2 + h * am2);

        s += h;

        if s > 0.049999 {
            s = 0.0;
            array_t.push(t);
            array_rm1.push(rm1);
            array_rm2.push(rm2);
            array_rt1.push(rt1);
            array_rt2.push(rt2);
            array_xncg.push(xnc / 32.2);
        }
    }

    Results {
        time: array_t,
        rm1: array_rm1,
        rm2: array_rm2,
        rt1: array_rt1,
        rt2: array_rt2,
        xncg: array_xncg,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c20l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.rm1.clone(),
        results.rm2.clone(),
        results.rt1.clone(),
        results.rt2.clone(),
        results.xncg.clone(),
    ])?;

    let plot_file = format!("{}/c20l1_trajectory.png", output_dir);
    let config = PlotConfig::new("Target Jink - Trajectories")
        .with_labels("Downrange (Ft)", "Altitude (Ft)");

    let series = vec![
        Series::new(results.rm1.clone(), results.rm2.clone())
            .with_color(plotters::prelude::BLUE)
            .with_label("Missile"),
        Series::new(results.rt1.clone(), results.rt2.clone())
            .with_color(plotters::prelude::RED)
            .with_label("Target"),
    ];

    line_plot(&plot_file, &config, &series).ok();

    // Compute miss distance
    let n = results.rm1.len();
    if n > 0 {
        let final_rtm = ((results.rt1[n-1] - results.rm1[n-1]).powi(2) +
                        (results.rt2[n-1] - results.rm2[n-1]).powi(2)).sqrt();
        println!("C20L1: Simulation finished");
        println!("  Data saved to: {}", data_file);
        println!("  Final miss distance: {:.2} ft", final_rtm);
    }

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c20l1_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }

    #[test]
    fn test_c20l1_small_miss() {
        let results = run();
        let n = results.rm1.len();
        if n > 0 {
            let final_rtm = ((results.rt1[n-1] - results.rm1[n-1]).powi(2) +
                            (results.rt2[n-1] - results.rm2[n-1]).powi(2)).sqrt();
            // With N'=3 and 200ft jink at 1 sec TGO, miss should be reasonable
            assert!(final_rtm < 300.0);
        }
    }
}
