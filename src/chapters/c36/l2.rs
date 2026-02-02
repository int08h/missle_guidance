//! Chapter 36, Lesson 2: Impact Angle Control 2D Engagement
//!
//! 2D engagement simulation with terminal impact angle constraint
//! against a non-maneuvering target.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub rt1k: Vec<f64>,
    pub rt2k: Vec<f64>,
    pub rm1k: Vec<f64>,
    pub rm2k: Vec<f64>,
    pub xncg: Vec<f64>,
    pub xlamdeg: Vec<f64>,
}

/// Run the C36L2 simulation
pub fn run() -> Results {
    let xntg: f64 = 0.0;
    let hedeg: f64 = 0.0;
    let xnp: f64 = 3.0;
    let rm1ic: f64 = 0.0;
    let rm2ic: f64 = 10000.0;
    let rt1ic: f64 = 30000.0;
    let rt2ic: f64 = 0.0;
    let vm: f64 = 3000.0;
    let vt: f64 = 0.0;
    let xnclimg: f64 = 9999999.0;
    let apn: i32 = 1;
    let xlamfdeg: f64 = -90.0;
    let mut h: f64;

    let xnclim = 32.2 * xnclimg;
    let xlamf = xlamfdeg / 57.3;
    let xnt = 32.2 * xntg;

    let mut rm1 = rm1ic;
    let mut rm2 = rm2ic;
    let rt1 = rt1ic;
    let rt2 = rt2ic;
    let mut beta: f64 = 0.0;
    let vt1 = -vt * beta.cos();
    let vt2 = vt * beta.sin();
    let he = hedeg / 57.3;

    let mut t: f64 = 0.0;
    let mut s: f64 = 0.0;
    let mut rtm1 = rt1 - rm1;
    let mut rtm2 = rt2 - rm2;
    let mut rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
    let mut xlam = rtm2.atan2(rtm1);
    let xlead = (vt * (beta + xlam).sin() / vm).asin();
    let thet = xlam + xlead;
    let mut vm1 = vm * (thet + he).cos();
    let mut vm2 = vm * (thet + he).sin();
    let vtm1 = vt1 - vm1;
    let vtm2 = vt2 - vm2;
    let mut vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;

    let mut array_t = Vec::new();
    let mut array_rt1k = Vec::new();
    let mut array_rt2k = Vec::new();
    let mut array_rm1k = Vec::new();
    let mut array_rm2k = Vec::new();
    let mut array_xncg = Vec::new();
    let mut array_xlamdeg = Vec::new();

    while vc >= 0.0 {
        if rtm < 1000.0 {
            h = 0.00001;
        } else {
            h = 0.0001;
        }

        let betaold = beta;
        let _rt1old = rt1;
        let _rt2old = rt2;
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
        xlam = rtm2.atan2(rtm1);
        let xlamd = (rtm1 * vtm2 - rtm2 * vtm1) / (rtm * rtm);
        let tgo = rtm / vc;

        let mut xnc = if apn == 0 {
            xnp * vc * xlamd
        } else {
            let xnt1 = xnt * beta.sin();
            let xnt2 = xnt * beta.cos();
            let xntplos = -xnt1 * xlam.sin() + xnt2 * xlam.cos();
            4.0 * vc * xlamd + xntplos + 2.0 * vc * (xlam - xlamf) / tgo
        };

        if xnc > xnclim {
            xnc = xnclim;
        }
        if xnc < -xnclim {
            xnc = -xnclim;
        }

        let am1 = -xnc * xlam.sin();
        let am2 = xnc * xlam.cos();
        let _vt1 = -vt * beta.cos();
        let _vt2 = vt * beta.sin();
        let betad = if vt == 0.0 { 0.0 } else { xnt / vt };

        // Euler step
        beta += h * betad;
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
        xlam = rtm2.atan2(rtm1);
        let xlamd = (rtm1 * vtm2 - rtm2 * vtm1) / (rtm * rtm);
        let tgo = rtm / vc;

        let mut xnc = if apn == 0 {
            xnp * vc * xlamd
        } else {
            let xnt1 = xnt * beta.sin();
            let xnt2 = xnt * beta.cos();
            let xntplos = -xnt1 * xlam.sin() + xnt2 * xlam.cos();
            4.0 * vc * xlamd + xntplos + 2.0 * vc * (xlam - xlamf) / tgo
        };

        if xnc > xnclim {
            xnc = xnclim;
        }
        if xnc < -xnclim {
            xnc = -xnclim;
        }

        let am1 = -xnc * xlam.sin();
        let am2 = xnc * xlam.cos();
        let betad = if vt == 0.0 { 0.0 } else { xnt / vt };

        // RK2 averaging
        beta = 0.5 * (betaold + beta + h * betad);
        rm1 = 0.5 * (rm1old + rm1 + h * vm1);
        rm2 = 0.5 * (rm2old + rm2 + h * vm2);
        vm1 = 0.5 * (vm1old + vm1 + h * am1);
        vm2 = 0.5 * (vm2old + vm2 + h * am2);

        s += h;
        if s >= 0.09999 {
            s = 0.0;
            let rt1k = rt1 / 1000.0;
            let rt2k = rt2 / 1000.0;
            let rm1k = rm1 / 1000.0;
            let rm2k = rm2 / 1000.0;
            let xlamdeg = xlam * 57.3;
            let xncg = xnc / 32.2;

            array_t.push(t);
            array_rt1k.push(rt1k);
            array_rt2k.push(rt2k);
            array_rm1k.push(rm1k);
            array_rm2k.push(rm2k);
            array_xncg.push(xncg);
            array_xlamdeg.push(xlamdeg);
        }
    }

    Results {
        time: array_t,
        rt1k: array_rt1k,
        rt2k: array_rt2k,
        rm1k: array_rm1k,
        rm2k: array_rm2k,
        xncg: array_xncg,
        xlamdeg: array_xlamdeg,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c36l2_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.rt1k.clone(),
        results.rt2k.clone(),
        results.rm1k.clone(),
        results.rm2k.clone(),
        results.xncg.clone(),
        results.xlamdeg.clone(),
    ])?;

    let plot_file = format!("{}/c36l2_traj.png", output_dir);
    let config = PlotConfig::new("Engagement Geometry")
        .with_labels("Downrange (Kft)", "Altitude (Kft)");

    let series = vec![
        Series::new(results.rt1k.clone(), results.rt2k.clone())
            .with_color(plotters::prelude::RED)
            .with_label("Target"),
        Series::new(results.rm1k.clone(), results.rm2k.clone())
            .with_color(plotters::prelude::BLUE)
            .with_label("Missile"),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C36L2: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c36l2_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
