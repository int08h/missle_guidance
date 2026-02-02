//! Chapter 32, Lesson 2: Ballistic Target Intercept with Predictive Guidance
//!
//! Simulates a missile intercepting a stationary target using predictive
//! guidance with atmospheric drag and gravity effects.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub rm1k: Vec<f64>,
    pub rm2k: Vec<f64>,
    pub rt1k: Vec<f64>,
    pub rt2k: Vec<f64>,
    pub xncg: Vec<f64>,
}

/// PREDICTP function - predict final downrange position
#[allow(clippy::too_many_arguments)]
fn predictp(tp: f64, rm1p: f64, rm2p: f64, vm1p: f64, vm2p: f64, xnc1f: f64, rt1p: f64, rt2p: f64, betah: f64) -> f64 {
    let mut rm1 = rm1p;
    let mut rm2 = rm2p;
    let mut vm1 = vm1p;
    let mut vm2 = vm2p;
    let xnc = xnc1f;
    let rt1 = rt1p;
    let rt2 = rt2p;
    let beta = betah;
    let mut _t = tp;

    let mut rtm1 = rt1 - rm1;
    let mut rtm2 = rt2 - rm2;
    let mut rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();

    while rm2 >= 0.0 {
        let h = if rtm < 1000.0 { 0.001 } else { 0.1 };

        let rm1old = rm1;
        let rm2old = rm2;
        let vm1old = vm1;
        let vm2old = vm2;

        // First derivative evaluation
        rtm1 = rt1 - rm1;
        rtm2 = rt2 - rm2;
        let _rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
        let xlam = rtm2.atan2(rtm1);
        let xne1 = -xnc * xlam.sin();
        let xne2 = xnc * xlam.cos();
        let rho = if rm2 < 30000.0 {
            0.002378 * (-rm2 / 30000.0).exp()
        } else {
            0.0034 * (-rm2 / 22000.0).exp()
        };
        let vm = (vm1 * vm1 + vm2 * vm2).sqrt();
        let q = 0.5 * rho * vm * vm;
        let gam = vm2.atan2(vm1);
        let drag = q * 32.2 / beta;
        let am1 = -drag * gam.cos() + xne1;
        let am2 = -32.2 - drag * gam.sin() + xne2;

        // Euler step
        rm1 += h * vm1;
        rm2 += h * vm2;
        vm1 += h * am1;
        vm2 += h * am2;
        _t += h;

        // Second derivative evaluation
        rtm1 = rt1 - rm1;
        rtm2 = rt2 - rm2;
        rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
        let xlam = rtm2.atan2(rtm1);
        let xne1 = -xnc * xlam.sin();
        let xne2 = xnc * xlam.cos();
        let rho = if rm2 < 30000.0 {
            0.002378 * (-rm2 / 30000.0).exp()
        } else {
            0.0034 * (-rm2 / 22000.0).exp()
        };
        let vm = (vm1 * vm1 + vm2 * vm2).sqrt();
        let q = 0.5 * rho * vm * vm;
        let gam = vm2.atan2(vm1);
        let drag = q * 32.2 / beta;
        let am1 = -drag * gam.cos() + xne1;
        let am2 = -32.2 - drag * gam.sin() + xne2;

        // RK2 averaging
        rm1 = 0.5 * (rm1old + rm1 + h * vm1);
        rm2 = 0.5 * (rm2old + rm2 + h * vm2);
        vm1 = 0.5 * (vm1old + vm1 + h * am1);
        vm2 = 0.5 * (vm2old + vm2 + h * am2);
    }

    rm1
}

/// Run the C32L2 simulation
pub fn run() -> Results {
    let vm_init: f64 = 3000.0;
    let beta: f64 = 1000.0;
    let betah: f64 = 1000.0;
    let xlimg: f64 = 5.0;
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
    let mut array_rm1k = Vec::new();
    let mut array_rm2k = Vec::new();
    let mut array_rt1k = Vec::new();
    let mut array_rt2k = Vec::new();
    let mut array_xncg = Vec::new();

    while !(t > 0.0 && rm2 <= 0.0) {
        let rtm1 = rt1 - rm1;
        let rtm2 = rt2 - rm2;
        let rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();

        let h = if rtm < 1000.0 { 0.001 } else { 0.1 };

        let rm1old = rm1;
        let rm2old = rm2;
        let vm1old = vm1;
        let vm2old = vm2;

        // First derivative evaluation
        let rtm1 = rt1 - rm1;
        let rtm2 = rt2 - rm2;
        let _rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
        let xlam = rtm2.atan2(rtm1);
        let xne1 = -xnc * xlam.sin();
        let xne2 = xnc * xlam.cos();
        let rho = if rm2 < 30000.0 {
            0.002378 * (-rm2 / 30000.0).exp()
        } else {
            0.0034 * (-rm2 / 22000.0).exp()
        };
        let vm = (vm1 * vm1 + vm2 * vm2).sqrt();
        let q = 0.5 * rho * vm * vm;
        let gam = vm2.atan2(vm1);
        let drag = q * 32.2 / beta;
        let am1 = -drag * gam.cos() + xne1;
        let am2 = -32.2 - drag * gam.sin() + xne2;

        // Euler step
        rm1 += h * vm1;
        rm2 += h * vm2;
        vm1 += h * am1;
        vm2 += h * am2;
        t += h;

        // Second derivative evaluation
        let rtm1 = rt1 - rm1;
        let rtm2 = rt2 - rm2;
        let _rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
        let xlam = rtm2.atan2(rtm1);
        let xne1 = -xnc * xlam.sin();
        let xne2 = xnc * xlam.cos();
        let rho = if rm2 < 30000.0 {
            0.002378 * (-rm2 / 30000.0).exp()
        } else {
            0.0034 * (-rm2 / 22000.0).exp()
        };
        let vm = (vm1 * vm1 + vm2 * vm2).sqrt();
        let q = 0.5 * rho * vm * vm;
        let gam = vm2.atan2(vm1);
        let drag = q * 32.2 / beta;
        let am1 = -drag * gam.cos() + xne1;
        let am2 = -32.2 - drag * gam.sin() + xne2;

        // RK2 averaging
        rm1 = 0.5 * (rm1old + rm1 + h * vm1);
        rm2 = 0.5 * (rm2old + rm2 + h * vm2);
        vm1 = 0.5 * (vm1old + vm1 + h * am1);
        vm2 = 0.5 * (vm2old + vm2 + h * am2);

        s += h;
        if s >= 0.99999 {
            s = 0.0;
            if t > 30.0 {
                let x1 = predictp(t, rm1, rm2, vm1, vm2, xnc, rt1, rt2, betah);
                let delx = rt1 - x1;
                let xncp = xnc + 1.0;
                let x2 = predictp(t, rm1, rm2, vm1, vm2, xncp, rt1, rt2, betah);
                let dxdnc = (x2 - x1) / (xncp - xnc);
                let delxnc = delx / dxdnc;
                xnc += delxnc;
                if xnc > xlim {
                    xnc = xlim;
                }
                if xnc < -xlim {
                    xnc = -xlim;
                }
            }

            let rm1k = rm1 / 1000.0;
            let rm2k = rm2 / 1000.0;
            let rt1k = rt1 / 1000.0;
            let rt2k = rt2 / 1000.0;
            let xncg = xnc / 32.2;

            array_t.push(t);
            array_rm1k.push(rm1k);
            array_rm2k.push(rm2k);
            array_rt1k.push(rt1k);
            array_rt2k.push(rt2k);
            array_xncg.push(xncg);
        }
    }

    Results {
        time: array_t,
        rm1k: array_rm1k,
        rm2k: array_rm2k,
        rt1k: array_rt1k,
        rt2k: array_rt2k,
        xncg: array_xncg,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c32l2_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.rm1k.clone(),
        results.rm2k.clone(),
        results.rt1k.clone(),
        results.rt2k.clone(),
        results.xncg.clone(),
    ])?;

    let plot_file = format!("{}/c32l2_traj.png", output_dir);
    let config = PlotConfig::new("Ballistic Target Intercept with Predictive Guidance")
        .with_labels("Downrange (Kft)", "Altitude (Kft)");

    let series = vec![
        Series::new(results.rm1k.clone(), results.rm2k.clone())
            .with_color(plotters::prelude::BLUE)
            .with_label("Missile"),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C32L2: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c32l2_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
