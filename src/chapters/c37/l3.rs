//! Chapter 37, Lesson 3: Dive Guidance with Polynomial Trajectory
//!
//! 2D engagement with polynomial guidance for terminal dive.
//!
//! MATLAB BUG: Line 44 has `while RM@ >= 0` which should be `while RM2 >= 0`.
//! The `@` is a typo for `2`. This implementation uses the correct condition.
//! Output: 8 columns [T, RT1K, RT2K, RM1K, RM2K, XNCG, EPSDEG, RTM]

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub rt1k: Vec<f64>,
    pub rt2k: Vec<f64>,
    pub rm1k: Vec<f64>,
    pub rm2k: Vec<f64>,
    pub xncg: Vec<f64>,
    pub epsdeg: Vec<f64>,
    pub rtm: Vec<f64>,
}

/// Run the C37L3 simulation
pub fn run() -> Results {
    let iguid: i32 = 2;
    let xn: f64 = 4.0;  // Note: MATLAB has both XN=3 and XN=4, last one wins
    let tfdes: f64 = 50.0;
    let gamdeg: f64 = 30.0;
    let gam_init = gamdeg / 57.3;
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
    let mut h: f64;

    let xnclim = 32.2 * xnclimg;
    let mut rm1 = rm1ic;
    let mut rm2 = rm2ic;
    let rt1 = rt1ic;
    let rt2 = rt2ic;
    let vt1: f64 = 0.0;
    let vt2: f64 = 0.0;

    let mut t: f64 = 0.0;
    let mut s: f64 = 0.0;
    let mut rtm1 = rt1 - rm1;
    let mut rtm2 = rt2 - rm2;
    let mut rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
    let mut xlam = rtm2.atan2(rtm1);
    let mut vm1 = vm_init * gam_init.cos();
    let mut vm2 = vm_init * gam_init.sin();
    let vtm1 = vt1 - vm1;
    let vtm2 = vt2 - vm2;
    let _vc_init = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;
    let mut vc: f64;
    let gamf = gamfdeg / 57.3;
    let biasdeg = (-gamfdeg * (xnp - 1.0) + xnp * xlam * 57.3 - gamdeg) / delt;
    let bias = biasdeg / 57.3;
    let mut x: f64 = 0.0;
    let mut gam = gam_init;

    let mut array_t = Vec::new();
    let mut array_rt1k = Vec::new();
    let mut array_rt2k = Vec::new();
    let mut array_rm1k = Vec::new();
    let mut array_rm2k = Vec::new();
    let mut array_xncg = Vec::new();
    let mut array_epsdeg = Vec::new();
    let mut array_rtm = Vec::new();

    // Loop condition: RM2 >= 0 (missile altitude above ground)
    while rm2 >= 0.0 {
        if rtm < 1000.0 {
            h = 0.000001;
        } else {
            h = 0.001;
        }

        let gamold = gam;
        let xlamold = xlam;
        let rm1old = rm1;
        let rm2old = rm2;
        let vm1old = vm1;
        let vm2old = vm2;
        let xold = x;

        // First derivative evaluation
        let vm = (vm1 * vm1 + vm2 * vm2).sqrt();
        rtm1 = rt1 - rm1;
        rtm2 = rt2 - rm2;
        rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
        let vtm1 = vt1 - vm1;
        let vtm2 = vt2 - vm2;
        vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;
        let xlamd = (rtm1 * vtm2 - rtm2 * vtm1) / (rtm * rtm);

        let (mut xnc, mut gamd) = if iguid == 0 {
            let tgo = rtm / vc;
            let xnc = 4.0 * vc * xlamd + 2.0 * vc * (xlam - gamf) / tgo;
            (xnc, xnc / vm)
        } else if iguid == 1 {
            if t < tbeg {
                let gamd = xnp * xlamd;
                (vm * gamd, gamd)
            } else if t < tend {
                let gamd = xnp * xlamd + bias;
                (vm * gamd, gamd)
            } else {
                let gamd = xnp * xlamd;
                (vm * gamd, gamd)
            }
        } else {
            // iguid == 2: polynomial guidance
            let eps = gam - xlam;
            let temp1 = (xn - 2.0 + 2.0 * eps.cos()) * vm * (tfdes - t) - xn * rtm;
            let x2 = 0.5 * (xn - 1.0) * temp1 / (tfdes - t).powi(2);
            let gamd = 2.0 * x2 / (vm * eps.sin()) + xlamd;
            (vm * gamd, gamd)
        };

        if xnc > xnclim {
            xnc = xnclim;
            gamd = xnclim / vm;
        }
        if xnc < -xnclim {
            xnc = -xnclim;
            gamd = -xnclim / vm;
        }

        let _rho = 0.002378 * (-rm2 / 22000.0).exp();
        let am1 = -xnc * gam.sin();
        let am2 = xnc * gam.cos();
        let xd = xnc * xnc;

        // Euler step
        gam += h * gamd;
        xlam += h * xlamd;
        rm1 += h * vm1;
        rm2 += h * vm2;
        vm1 += h * am1;
        vm2 += h * am2;
        x += h * xd;
        t += h;

        // Second derivative for RK2
        let vm = (vm1 * vm1 + vm2 * vm2).sqrt();
        rtm1 = rt1 - rm1;
        rtm2 = rt2 - rm2;
        rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
        let vtm1 = vt1 - vm1;
        let vtm2 = vt2 - vm2;
        vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;
        let xlamd = (rtm1 * vtm2 - rtm2 * vtm1) / (rtm * rtm);

        let (mut xnc, mut gamd) = if iguid == 0 {
            let tgo = rtm / vc;
            let xnc = 4.0 * vc * xlamd + 2.0 * vc * (xlam - gamf) / tgo;
            (xnc, xnc / vm)
        } else if iguid == 1 {
            if t < tbeg {
                let gamd = xnp * xlamd;
                (vm * gamd, gamd)
            } else if t < tend {
                let gamd = xnp * xlamd + bias;
                (vm * gamd, gamd)
            } else {
                let gamd = xnp * xlamd;
                (vm * gamd, gamd)
            }
        } else {
            let eps = gam - xlam;
            let temp1 = (xn - 2.0 + 2.0 * eps.cos()) * vm * (tfdes - t) - xn * rtm;
            let x2 = 0.5 * (xn - 1.0) * temp1 / (tfdes - t).powi(2);
            let gamd = 2.0 * x2 / (vm * eps.sin()) + xlamd;
            (vm * gamd, gamd)
        };

        if xnc > xnclim {
            xnc = xnclim;
            gamd = xnclim / vm;
        }
        if xnc < -xnclim {
            xnc = -xnclim;
            gamd = -xnclim / vm;
        }

        let am1 = -xnc * gam.sin();
        let am2 = xnc * gam.cos();
        let xd = xnc * xnc;

        // RK2 averaging
        gam = 0.5 * (gamold + gam + h * gamd);
        xlam = 0.5 * (xlamold + xlam + h * xlamd);
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
            let _gamdeg = gam * 57.3;
            let _vmm = vm / 3.28;
            let _xm = x / (3.28 * 3.28);
            let eps = gam - xlam;
            let epsdeg = eps * 57.3;

            array_t.push(t);
            array_rt1k.push(rt1k);
            array_rt2k.push(rt2k);
            array_rm1k.push(rm1k);
            array_rm2k.push(rm2k);
            array_xncg.push(xncg);
            array_epsdeg.push(epsdeg);
            array_rtm.push(rtm);
        }
    }

    Results {
        time: array_t,
        rt1k: array_rt1k,
        rt2k: array_rt2k,
        rm1k: array_rm1k,
        rm2k: array_rm2k,
        xncg: array_xncg,
        epsdeg: array_epsdeg,
        rtm: array_rtm,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c37l3_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.rt1k.clone(),
        results.rt2k.clone(),
        results.rm1k.clone(),
        results.rm2k.clone(),
        results.xncg.clone(),
        results.epsdeg.clone(),
        results.rtm.clone(),
    ])?;

    let plot_file = format!("{}/c37l3_traj.png", output_dir);
    let config = PlotConfig::new("Engagement Geometry")
        .with_labels("Downrange (Km)", "Altitude (Km)");

    let series = vec![
        Series::new(results.rt1k.clone(), results.rt2k.clone())
            .with_color(plotters::prelude::RED)
            .with_label("Target"),
        Series::new(results.rm1k.clone(), results.rm2k.clone())
            .with_color(plotters::prelude::BLUE)
            .with_label("Missile"),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C37L3: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c37l3_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
