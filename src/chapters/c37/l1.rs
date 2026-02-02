//! Chapter 37, Lesson 1: Dive Guidance
//!
//! Terminal dive guidance with flight path angle control.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub rt1k: Vec<f64>,
    pub rt2k: Vec<f64>,
    pub rm1k: Vec<f64>,
    pub rm2k: Vec<f64>,
    pub xncg: Vec<f64>,
    pub gamdeg: Vec<f64>,
}

/// Run the C37L1 simulation
pub fn run() -> Results {
    let iguid: i32 = 0;
    let ichoice: i32 = 0;
    let gamdeg_init: f64 = 30.0;
    let gamic = gamdeg_init / 57.3;
    let delt: f64 = 30.0;
    let tbeg: f64 = 0.0;
    let tend = tbeg + delt;
    let xnp: f64 = 3.0;
    let rm1ic: f64 = 0.0;
    let rm2ic: f64 = 0.0;
    let rt1ic: f64 = 10000.0 * 3.28;
    let rt2ic: f64 = 0.0;
    let vm_init: f64 = 250.0 * 3.28;
    let xnclimg: f64 = 10.0;
    let gamfdeg: f64 = -90.0;

    let xnclim = 32.2 * xnclimg;
    let mut rm1 = rm1ic;
    let mut rm2 = rm2ic;
    let rt1 = rt1ic;
    let rt2 = rt2ic;
    let vt1: f64 = 0.0;
    let vt2: f64 = 0.0;

    let mut t: f64 = 0.0;
    let mut s: f64 = 0.0;
    let mut h: f64;

    let mut rtm1 = rt1 - rm1;
    let mut rtm2 = rt2 - rm2;
    let mut rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
    let mut xlam = rtm2.atan2(rtm1);
    let mut vm1 = vm_init * gamic.cos();
    let mut vm2 = vm_init * gamic.sin();
    let vtm1 = vt1 - vm1;
    let vtm2 = vt2 - vm2;
    let mut vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;
    let gamf = gamfdeg / 57.3;
    let biasdeg = (-gamfdeg * (xnp - 1.0) + xnp * xlam * 57.3 - gamdeg_init) / delt;
    let bias = biasdeg / 57.3;
    let mut x: f64 = 0.0;

    let mut array_t = Vec::new();
    let mut array_rt1k = Vec::new();
    let mut array_rt2k = Vec::new();
    let mut array_rm1k = Vec::new();
    let mut array_rm2k = Vec::new();
    let mut array_xncg = Vec::new();
    let mut array_gamdeg = Vec::new();

    while vc >= 0.0 {
        if rtm < 1000.0 {
            h = 0.00001;
        } else {
            h = 0.0001;
        }

        let rm1old = rm1;
        let rm2old = rm2;
        let vm1old = vm1;
        let vm2old = vm2;
        let xold = x;

        // First derivative evaluation
        let gam = vm2.atan2(vm1);
        let vm = (vm1 * vm1 + vm2 * vm2).sqrt();
        rtm1 = rt1 - rm1;
        rtm2 = rt2 - rm2;
        rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
        let vtm1 = vt1 - vm1;
        let vtm2 = vt2 - vm2;
        vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;
        xlam = rtm2.atan2(rtm1);
        let xlamd = (rtm1 * vtm2 - rtm2 * vtm1) / (rtm * rtm);

        let mut xnc = if iguid == 0 {
            let tgo = rtm / vc;
            4.0 * vc * xlamd + 2.0 * vc * (xlam - gamf) / tgo
        } else if t < tbeg {
            let gamd = xnp * xlamd;
            vm * gamd
        } else if t < tend {
            let gamd = xnp * xlamd + bias;
            vm * gamd
        } else {
            let gamd = xnp * xlamd;
            vm * gamd
        };

        if xnc > xnclim {
            xnc = xnclim;
        }
        if xnc < -xnclim {
            xnc = -xnclim;
        }

        let (am1, am2) = if ichoice == 0 {
            (-xnc * xlam.sin(), xnc * xlam.cos())
        } else {
            (-xnc * gam.sin(), xnc * gam.cos())
        };
        let xd = xnc * xnc;

        // Euler step
        rm1 += h * vm1;
        rm2 += h * vm2;
        vm1 += h * am1;
        vm2 += h * am2;
        x += h * xd;
        t += h;

        // Second derivative for RK2
        let gam = vm2.atan2(vm1);
        let vm = (vm1 * vm1 + vm2 * vm2).sqrt();
        rtm1 = rt1 - rm1;
        rtm2 = rt2 - rm2;
        rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
        let vtm1 = vt1 - vm1;
        let vtm2 = vt2 - vm2;
        vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;
        xlam = rtm2.atan2(rtm1);
        let xlamd = (rtm1 * vtm2 - rtm2 * vtm1) / (rtm * rtm);

        let mut xnc = if iguid == 0 {
            let tgo = rtm / vc;
            4.0 * vc * xlamd + 2.0 * vc * (xlam - gamf) / tgo
        } else if t < tbeg {
            let gamd = xnp * xlamd;
            vm * gamd
        } else if t < tend {
            let gamd = xnp * xlamd + bias;
            vm * gamd
        } else {
            let gamd = xnp * xlamd;
            vm * gamd
        };

        if xnc > xnclim {
            xnc = xnclim;
        }
        if xnc < -xnclim {
            xnc = -xnclim;
        }

        let (am1, am2) = if ichoice == 0 {
            (-xnc * xlam.sin(), xnc * xlam.cos())
        } else {
            (-xnc * gam.sin(), xnc * gam.cos())
        };
        let xd = xnc * xnc;

        // RK2 averaging
        rm1 = 0.5 * (rm1old + rm1 + h * vm1);
        rm2 = 0.5 * (rm2old + rm2 + h * vm2);
        vm1 = 0.5 * (vm1old + vm1 + h * am1);
        vm2 = 0.5 * (vm2old + vm2 + h * am2);
        x = 0.5 * (xold + x + h * xd);

        s += h;
        if s >= 0.09999 {
            s = 0.0;
            let rt1k = rt1 / 3280.0;
            let rt2k = rt2 / 3280.0;
            let rm1k = rm1 / 3280.0;
            let rm2k = rm2 / 3280.0;
            let xncg = xnc / 32.2;
            let gamdeg = gam * 57.3;
            array_t.push(t);
            array_rt1k.push(rt1k);
            array_rt2k.push(rt2k);
            array_rm1k.push(rm1k);
            array_rm2k.push(rm2k);
            array_xncg.push(xncg);
            array_gamdeg.push(gamdeg);
        }
    }

    Results {
        time: array_t,
        rt1k: array_rt1k,
        rt2k: array_rt2k,
        rm1k: array_rm1k,
        rm2k: array_rm2k,
        xncg: array_xncg,
        gamdeg: array_gamdeg,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c37l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.rt1k.clone(),
        results.rt2k.clone(),
        results.rm1k.clone(),
        results.rm2k.clone(),
        results.xncg.clone(),
        results.gamdeg.clone(),
    ])?;

    let plot_file = format!("{}/c37l1_traj.png", output_dir);
    let config = PlotConfig::new("Dive Guidance - Trajectory")
        .with_labels("Downrange (Km)", "Altitude (Km)");

    let series = vec![
        Series::new(results.rm1k.clone(), results.rm2k.clone())
            .with_color(plotters::prelude::BLUE)
            .with_label("Missile"),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C37L1: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c37l1_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
