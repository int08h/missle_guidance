//! Chapter 11, Lesson 2: 2D Engagement with Reentry Target
//!
//! Two-dimensional engagement simulation with ballistic reentry target.
//!
//! MATLAB NOTE: The MATLAB code calls `initialpz` which is defined in INITIALPZ.m.
//! This function is implemented inline in Rust below.
//! Output: 8 columns [T, RT1K, RT2K, RM1K, RM2K, ATG, ATPLOSG, XNCG]

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub rt1k: Vec<f64>,
    pub rt2k: Vec<f64>,
    pub rm1k: Vec<f64>,
    pub rm2k: Vec<f64>,
    pub atg: Vec<f64>,
    pub atplosg: Vec<f64>,
    pub xncg: Vec<f64>,
}

/// Compute initial target position when it reaches desired altitude
fn initialpz(
    rt2des: f64,
    rt1ic: f64,
    rt2ic: f64,
    vt1ic: f64,
    vt2ic: f64,
    beta: f64,
) -> (f64, f64, f64) {
    let mut rt1 = rt1ic;
    let mut rt2 = rt2ic;
    let mut vt1 = vt1ic;
    let mut vt2 = vt2ic;
    let mut t: f64 = 0.0;
    let h: f64 = 0.01;

    while rt2 > rt2des {
        let rt1_old = rt1;
        let rt2_old = rt2;
        let vt1_old = vt1;
        let vt2_old = vt2;

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

        // Second derivative (at new state)
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
        rt1 = 0.5 * (rt1_old + rt1 + h * vt1);
        rt2 = 0.5 * (rt2_old + rt2 + h * vt2);
        vt1 = 0.5 * (vt1_old + vt1 + h * at1);
        vt2 = 0.5 * (vt2_old + vt2 + h * at2);
    }

    (rt1, rt2, t)
}

/// Run the C11L2 simulation
pub fn run() -> Results {
    let apn: i32 = 0;
    let xnp: f64 = 3.0;
    let mut rt1: f64 = 0.0;
    let mut rt2: f64 = 200000.0;
    let mut rm1: f64 = 170000.0;
    let mut rm2: f64 = 0.0;
    let vt: f64 = 6000.0;
    let rt2des: f64 = 50000.0;
    let gamtdeg: f64 = 45.0;
    let beta: f64 = 500.0;
    let betest: f64 = 500.0;
    let xnclimg: f64 = 7.0;
    let xnclim: f64 = xnclimg * 32.2;

    let mut vt1: f64 = vt * (gamtdeg / 57.3).cos();
    let mut vt2: f64 = -vt * (gamtdeg / 57.3).sin();

    // Compute initial conditions using target trajectory prediction
    let (rt1f, rt2f, tfdes) = initialpz(rt2des, rt1, rt2, vt1, vt2, betest);
    let rtm1f = rt1f - rm1;
    let rtm2f = rt2f - rm2;
    let gammdeg = 57.3 * rtm2f.atan2(rtm1f);
    let rtmf = (rtm1f * rtm1f + rtm2f * rtm2f).sqrt();
    let vm = rtmf / tfdes;
    let mut vm1 = vm * (gammdeg / 57.3).cos();
    let mut vm2 = vm * (gammdeg / 57.3).sin();

    let mut rtm1 = rt1 - rm1;
    let mut rtm2 = rt2 - rm2;
    let mut rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
    let vtm1 = vt1 - vm1;
    let vtm2 = vt2 - vm2;
    let mut vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;

    let mut t: f64 = 0.0;
    let mut h: f64;
    let mut s: f64 = 0.0;
    let mut xnc: f64;

    let mut array_t = Vec::new();
    let mut array_rt1k = Vec::new();
    let mut array_rt2k = Vec::new();
    let mut array_rm1k = Vec::new();
    let mut array_rm2k = Vec::new();
    let mut array_atg = Vec::new();
    let mut array_atplosg = Vec::new();
    let mut array_xncg = Vec::new();

    let mut at1: f64;
    let mut at2: f64;
    let mut am1: f64;
    let mut am2: f64;
    let mut atplos: f64;

    while vc >= 0.0 {
        if rtm < 1000.0 {
            h = 0.0002;
        } else {
            h = 0.01;
        }

        let rt1_old = rt1;
        let rt2_old = rt2;
        let vt1_old = vt1;
        let vt2_old = vt2;
        let rm1_old = rm1;
        let rm2_old = rm2;
        let vm1_old = vm1;
        let vm2_old = vm2;

        // First derivative evaluation
        let rho = if rt2 <= 30000.0 {
            0.002378 * (-rt2 / 30000.0).exp()
        } else {
            0.0034 * (-rt2 / 22000.0).exp()
        };
        let vt_mag = (vt1 * vt1 + vt2 * vt2).sqrt();
        let q = 0.5 * rho * vt_mag * vt_mag;
        let gamt = (-vt2).atan2(vt1);
        at1 = -32.2 * q * gamt.cos() / beta;
        at2 = -32.2 + 32.2 * q * gamt.sin() / beta;

        rtm1 = rt1 - rm1;
        rtm2 = rt2 - rm2;
        rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
        let vtm1 = vt1 - vm1;
        let vtm2 = vt2 - vm2;
        vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;
        let xlam = rtm2.atan2(rtm1);
        let xlamd = (rtm1 * vtm2 - rtm2 * vtm1) / (rtm * rtm);
        atplos = -at1 * xlam.sin() + at2 * xlam.cos();

        xnc = if apn == 1 {
            xnp * vc * xlamd + 0.5 * xnp * atplos
        } else {
            xnp * vc * xlamd
        };

        if xnc > xnclim {
            xnc = xnclim;
        }
        if xnc < -xnclim {
            xnc = -xnclim;
        }

        am1 = -xnc * xlam.sin();
        am2 = xnc * xlam.cos();

        // Euler step
        rt1 += h * vt1;
        rt2 += h * vt2;
        vt1 += h * at1;
        vt2 += h * at2;
        rm1 += h * vm1;
        rm2 += h * vm2;
        vm1 += h * am1;
        vm2 += h * am2;
        t += h;

        // Second derivative evaluation
        let rho = if rt2 <= 30000.0 {
            0.002378 * (-rt2 / 30000.0).exp()
        } else {
            0.0034 * (-rt2 / 22000.0).exp()
        };
        let vt_mag = (vt1 * vt1 + vt2 * vt2).sqrt();
        let q = 0.5 * rho * vt_mag * vt_mag;
        let gamt = (-vt2).atan2(vt1);
        at1 = -32.2 * q * gamt.cos() / beta;
        at2 = -32.2 + 32.2 * q * gamt.sin() / beta;

        rtm1 = rt1 - rm1;
        rtm2 = rt2 - rm2;
        rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
        let vtm1 = vt1 - vm1;
        let vtm2 = vt2 - vm2;
        vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;
        let xlam = rtm2.atan2(rtm1);
        let xlamd = (rtm1 * vtm2 - rtm2 * vtm1) / (rtm * rtm);
        atplos = -at1 * xlam.sin() + at2 * xlam.cos();

        xnc = if apn == 1 {
            xnp * vc * xlamd + 0.5 * xnp * atplos
        } else {
            xnp * vc * xlamd
        };

        if xnc > xnclim {
            xnc = xnclim;
        }
        if xnc < -xnclim {
            xnc = -xnclim;
        }

        am1 = -xnc * xlam.sin();
        am2 = xnc * xlam.cos();

        // RK2 averaging
        rt1 = 0.5 * (rt1_old + rt1 + h * vt1);
        rt2 = 0.5 * (rt2_old + rt2 + h * vt2);
        vt1 = 0.5 * (vt1_old + vt1 + h * at1);
        vt2 = 0.5 * (vt2_old + vt2 + h * at2);
        rm1 = 0.5 * (rm1_old + rm1 + h * vm1);
        rm2 = 0.5 * (rm2_old + rm2 + h * vm2);
        vm1 = 0.5 * (vm1_old + vm1 + h * am1);
        vm2 = 0.5 * (vm2_old + vm2 + h * am2);

        s += h;
        if s >= 0.09999 {
            s = 0.0;
            let atg = (at1 * at1 + at2 * at2).sqrt() / 32.2;
            let rt1k = rt1 / 1000.0;
            let rt2k = rt2 / 1000.0;
            let rm1k = rm1 / 1000.0;
            let rm2k = rm2 / 1000.0;
            let xncg = xnc / 32.2;
            let atplosg = atplos / 32.2;

            array_t.push(t);
            array_rt1k.push(rt1k);
            array_rt2k.push(rt2k);
            array_rm1k.push(rm1k);
            array_rm2k.push(rm2k);
            array_atg.push(atg);
            array_atplosg.push(atplosg);
            array_xncg.push(xncg);
        }
    }

    println!("RTM = {:.6}", rtm);

    Results {
        time: array_t,
        rt1k: array_rt1k,
        rt2k: array_rt2k,
        rm1k: array_rm1k,
        rm2k: array_rm2k,
        atg: array_atg,
        atplosg: array_atplosg,
        xncg: array_xncg,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c11l2_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.rt1k.clone(),
        results.rt2k.clone(),
        results.rm1k.clone(),
        results.rm2k.clone(),
        results.atg.clone(),
        results.atplosg.clone(),
        results.xncg.clone(),
    ])?;

    // Trajectory plot
    let plot_file = format!("{}/c11l2_trajectory.png", output_dir);
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
    let plot_file2 = format!("{}/c11l2_acceleration.png", output_dir);
    let config2 = PlotConfig::new("Accelerations")
        .with_labels("Time (Sec)", "Acceleration (G)");
    let series2 = vec![
        Series::new(results.time.clone(), results.atg.clone())
            .with_label("Target")
            .with_color(plotters::prelude::BLUE),
        Series::new(results.time.clone(), results.atplosg.clone())
            .with_label("ATPLOS")
            .with_color(plotters::prelude::GREEN),
        Series::new(results.time.clone(), results.xncg.clone())
            .with_label("Missile Cmd")
            .with_color(plotters::prelude::RED),
    ];
    line_plot(&plot_file2, &config2, &series2).ok();

    println!("C11L2: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c11l2_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
