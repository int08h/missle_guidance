//! Chapter 32, Lesson 1: Ballistic Target Intercept
//!
//! Simulates a missile intercepting a stationary target with atmospheric
//! drag and gravity effects.

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

/// Run the C32L1 simulation
pub fn run() -> Results {
    let ts: f64 = 0.1;
    let apn: i32 = 0;
    let xlimg: f64 = 5.0;
    let vm_init: f64 = 3000.0;
    let beta: f64 = 1000.0;
    let gamdeg: f64 = 45.0;

    let mut vm1 = vm_init * (gamdeg / 57.3).cos();
    let mut vm2 = vm_init * (gamdeg / 57.3).sin();
    let mut rm1: f64 = 0.0;
    let mut rm2: f64 = 0.0;
    let mut xnc: f64 = 0.0;
    let rt1: f64 = 60000.0;
    let rt2: f64 = 0.0;

    let mut t: f64 = 0.0;
    let mut s: f64 = 0.0;
    let xlim = xlimg * 32.2;

    let mut array_t = Vec::new();
    let mut array_rt1k = Vec::new();
    let mut array_rt2k = Vec::new();
    let mut array_rm1k = Vec::new();
    let mut array_rm2k = Vec::new();
    let mut array_xncg = Vec::new();

    while !(t > 0.0 && rm2 <= 0.0) {
        let h = if (rt1 - rm1).hypot(rt2 - rm2) < 1000.0 {
            0.0001
        } else {
            0.01
        };

        let rm1old = rm1;
        let rm2old = rm2;
        let vm1old = vm1;
        let vm2old = vm2;

        // First derivative evaluation
        let rtm1 = rt1 - rm1;
        let rtm2 = rt2 - rm2;
        let rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
        let vtm1 = -vm1;
        let vtm2 = -vm2;
        let _vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;
        let xlam = rtm2.atan2(rtm1);
        let _xlamd = (rtm1 * vtm2 - rtm2 * vtm1) / (rtm * rtm);

        let rho = if rm2 < 30000.0 {
            0.002378 * (-rm2 / 30000.0).exp()
        } else {
            0.0034 * (-rm2 / 22000.0).exp()
        };

        let vm = (vm1 * vm1 + vm2 * vm2).sqrt();
        let q = 0.5 * rho * vm * vm;
        let gam = vm2.atan2(vm1);
        let drag = q * 32.2 / beta;
        let xne1 = -xnc * xlam.sin();
        let xne2 = xnc * xlam.cos();
        let am1 = -drag * gam.cos() + xne1;
        let am2 = -32.2 - drag * gam.sin() + xne2;

        // Euler step
        rm1 += h * vm1;
        rm2 += h * vm2;
        vm1 += h * am1;
        vm2 += h * am2;
        t += h;

        // Second derivative for RK2
        let rtm1 = rt1 - rm1;
        let rtm2 = rt2 - rm2;
        let rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
        let vtm1 = -vm1;
        let vtm2 = -vm2;
        let vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;
        let xlam = rtm2.atan2(rtm1);
        let xlamd = (rtm1 * vtm2 - rtm2 * vtm1) / (rtm * rtm);

        let rho = if rm2 < 30000.0 {
            0.002378 * (-rm2 / 30000.0).exp()
        } else {
            0.0034 * (-rm2 / 22000.0).exp()
        };

        let vm = (vm1 * vm1 + vm2 * vm2).sqrt();
        let q = 0.5 * rho * vm * vm;
        let gam = vm2.atan2(vm1);
        let drag = q * 32.2 / beta;
        let xne1 = -xnc * xlam.sin();
        let xne2 = xnc * xlam.cos();
        let am1 = -drag * gam.cos() + xne1;
        let am2 = -32.2 - drag * gam.sin() + xne2;

        // RK2 averaging
        rm1 = 0.5 * (rm1old + rm1 + h * vm1);
        rm2 = 0.5 * (rm2old + rm2 + h * vm2);
        vm1 = 0.5 * (vm1old + vm1 + h * am1);
        vm2 = 0.5 * (vm2old + vm2 + h * am2);

        s += h;
        if s >= ts - 0.00001 {
            s = 0.0;
            let rm1k = rm1 / 1000.0;
            let rm2k = rm2 / 1000.0;
            let rt1k = rt1 / 1000.0;
            let rt2k = rt2 / 1000.0;

            let _drag1 = -drag * gam.cos();
            let drag2 = -drag * gam.sin() - 32.2;
            let dragplos = -(-drag * gam.cos()) * xlam.sin() + drag2 * xlam.cos();
            let atplos = 0.0;

            if t > 30.0 {
                if apn == 0 {
                    xnc = 3.0 * vc * xlamd;
                } else {
                    xnc = 3.0 * vc * xlamd + 1.5 * (atplos - dragplos);
                }
            } else {
                xnc = 0.0;
            }

            if xnc > xlim {
                xnc = xlim;
            }
            if xnc < -xlim {
                xnc = -xlim;
            }

            let xncg = xnc / 32.2;
            array_t.push(t);
            array_rt1k.push(rt1k);
            array_rt2k.push(rt2k);
            array_rm1k.push(rm1k);
            array_rm2k.push(rm2k);
            array_xncg.push(xncg);
        }
    }

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

    let data_file = format!("{}/c32l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.rt1k.clone(),
        results.rt2k.clone(),
        results.rm1k.clone(),
        results.rm2k.clone(),
        results.xncg.clone(),
    ])?;

    let plot_file = format!("{}/c32l1_traj.png", output_dir);
    let config = PlotConfig::new("Ballistic Target Intercept")
        .with_labels("Downrange (Kft)", "Altitude (Kft)");

    let series = vec![
        Series::new(results.rm1k.clone(), results.rm2k.clone())
            .with_color(plotters::prelude::BLUE)
            .with_label("Missile"),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C32L1: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c32l1_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
