//! Chapter 14, Lesson 1: 2D Engagement with Noise
//!
//! Full 2D missile-target engagement with proportional navigation,
//! fading memory filter for LOS rate estimation, and measurement noise.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;
use rand::prelude::*;
use rand::SeedableRng;
use rand_distr::StandardNormal;

pub struct Results {
    pub time: Vec<f64>,
    pub xlamd: Vec<f64>,
    pub xlamdh: Vec<f64>,
    pub xncg: Vec<f64>,
    pub rtm: Vec<f64>,
}

/// Run the C14L1 simulation
pub fn run() -> Results {
    run_with_seed(None)
}

/// Run with optional seed for reproducibility
pub fn run_with_seed(seed: Option<u64>) -> Results {
    let vm: f64 = 3000.0;
    let vt: f64 = 1000.0;
    let xnt: f64 = 0.0;
    let rm1ic: f64 = 0.001;
    let rm2ic: f64 = 10000.0;
    let rt1ic: f64 = 40000.0;
    let rt2ic: f64 = 10000.0;
    let hedeg: f64 = 0.0;
    let xnp: f64 = 3.0;
    let beta: f64 = 0.8;
    let ts: f64 = 0.1;
    let signoise: f64 = 0.001;
    let noise: bool = true;

    let mut rng: Box<dyn RngCore> = match seed {
        Some(s) => Box::new(rand::rngs::StdRng::seed_from_u64(s)),
        None => Box::new(rand::rngs::StdRng::from_entropy()),
    };

    let mut rt1 = rt1ic;
    let mut rt2 = rt2ic;
    let mut rm1 = rm1ic;
    let mut rm2 = rm2ic;
    let mut betat: f64 = 0.0;

    let mut vt1 = -vt * betat.cos();
    let mut vt2 = vt * betat.sin();
    let he = hedeg / 57.3;

    let gfilter = 1.0 - beta * beta;
    let hfilter = (1.0 - beta) * (1.0 - beta);

    let mut xlamh: f64 = 0.0;
    let mut xlamdh: f64 = 0.0;
    let mut xnc: f64;

    let mut t: f64 = 0.0;
    let mut s: f64 = 0.0;

    // Initial geometry
    let mut rtm1 = rt1 - rm1;
    let mut rtm2 = rt2 - rm2;
    let mut rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
    let xlam = rtm2.atan2(rtm1);
    let xlead = (vt * (betat + xlam).sin() / vm).asin();
    let thet = xlam + xlead;

    let mut vm1 = vm * (thet + he).cos();
    let mut vm2 = vm * (thet + he).sin();

    let vtm1 = vt1 - vm1;
    let vtm2 = vt2 - vm2;
    let mut vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;

    let mut array_t = Vec::new();
    let mut array_xlamd = Vec::new();
    let mut array_xlamdh = Vec::new();
    let mut array_xncg = Vec::new();
    let mut array_rtm = Vec::new();

    while vc >= 0.0 {
        let h = if rtm < 1000.0 { 0.0002 } else { 0.01 };

        let betatold = betat;
        let rt1old = rt1;
        let rt2old = rt2;
        let rm1old = rm1;
        let rm2old = rm2;
        let vm1old = vm1;
        let vm2old = vm2;

        // First derivative evaluation
        let thett = rt2.atan2(rt1);
        let thetm = rm2.atan2(rm1);
        let _rt = (rt1 * rt1 + rt2 * rt2).sqrt();
        let _rm = (rm1 * rm1 + rm2 * rm2).sqrt();

        rtm1 = rt1 - rm1;
        rtm2 = rt2 - rm2;
        rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();

        let vtm1 = vt1 - vm1;
        let vtm2 = vt2 - vm2;
        vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;

        let xlam = (rtm2 / rtm1).atan();
        let xlamd = (rtm1 * vtm2 - rtm2 * vtm1) / (rtm * rtm);

        xnc = xnp * vc * xlamdh;
        let am1 = -xnc * xlam.sin();
        let am2 = xnc * xlam.cos();

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

        // Second derivative for RK2
        let _thett_new = rt2.atan2(rt1);
        let _thetm_new = rm2.atan2(rm1);
        let _rt_new = (rt1 * rt1 + rt2 * rt2).sqrt();
        let _rm_new = (rm1 * rm1 + rm2 * rm2).sqrt();

        rtm1 = rt1 - rm1;
        rtm2 = rt2 - rm2;
        rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();

        let vtm1 = vt1 - vm1;
        let vtm2 = vt2 - vm2;
        vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;

        let xlam = (rtm2 / rtm1).atan();

        xnc = xnp * vc * xlamdh;
        let am1 = -xnc * xlam.sin();
        let am2 = xnc * xlam.cos();

        vt1 = -vt * betat.cos();
        vt2 = vt * betat.sin();
        let betatd = xnt / vt;

        // RK2 averaging
        betat = 0.5 * (betatold + betat + h * betatd);
        rt1 = 0.5 * (rt1old + rt1 + h * vt1);
        rt2 = 0.5 * (rt2old + rt2 + h * vt2);
        rm1 = 0.5 * (rm1old + rm1 + h * vm1);
        rm2 = 0.5 * (rm2old + rm2 + h * vm2);
        vm1 = 0.5 * (vm1old + vm1 + h * am1);
        vm2 = 0.5 * (vm2old + vm2 + h * am2);

        s += h;

        if s >= ts - 1e-5 {
            s = 0.0;

            let thettnoise = if noise {
                signoise * rng.sample::<f64, _>(StandardNormal)
            } else {
                0.0
            };

            let thettm = thett + thettnoise;
            let thetmm = thetm;

            let rt = (rt1 * rt1 + rt2 * rt2).sqrt();
            let rm = (rm1 * rm1 + rm2 * rm2).sqrt();

            let rt1m = rt * thettm.cos();
            let rt2m = rt * thettm.sin();
            let rm1m = rm * thetmm.cos();
            let rm2m = rm * thetmm.sin();

            let xlamm = (rt2m - rm2m).atan2(rt1m - rm1m);

            let res = xlamm - (xlamh + ts * xlamdh);
            xlamh = gfilter * res + xlamh + ts * xlamdh;
            xlamdh += hfilter * res / ts;

            xnc = xnp * vc * xlamdh;

            array_t.push(t);
            array_xlamd.push(xlamd);
            array_xlamdh.push(xlamdh);
            array_xncg.push(xnc / 32.2);
            array_rtm.push(rtm);
        }
    }

    Results {
        time: array_t,
        xlamd: array_xlamd,
        xlamdh: array_xlamdh,
        xncg: array_xncg,
        rtm: array_rtm,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c14l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.xlamd.clone(),
        results.xlamdh.clone(),
        results.xncg.clone(),
        results.rtm.clone(),
    ])?;

    let plot_file = format!("{}/c14l1_plot.png", output_dir);
    let config = PlotConfig::new("2D Engagement - LOS Rate Estimate")
        .with_labels("Time (S)", "Line of Sight Rate Estimate (Rad/S)")
        .with_y_range(-0.02, 0.02);

    let series = vec![
        Series::new(results.time.clone(), results.xlamdh.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C14L1: Simulation finished");
    println!("  Data saved to: {}", data_file);
    println!("  Plot saved to: {}", plot_file);
    println!("  Final miss distance: {:.2} ft", results.rtm.last().unwrap_or(&0.0));

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c14l1_runs() {
        let results = run_with_seed(Some(12345));
        assert!(!results.time.is_empty());
    }

    #[test]
    fn test_c14l1_intercept() {
        let results = run_with_seed(Some(12345));
        // Should achieve small miss distance
        let min_rtm = results.rtm.iter().cloned().fold(f64::MAX, f64::min);
        assert!(min_rtm < 100.0);  // Within 100 ft miss
    }
}
