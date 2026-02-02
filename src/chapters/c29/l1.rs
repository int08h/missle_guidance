//! Chapter 29, Lesson 1: 2D Engagement with Weaving Target
//!
//! Simulates miss distance for various initial ranges with a sinusoidally
//! weaving target.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub miss: Vec<f64>,
}

/// Run the C29L1 simulation
pub fn run() -> Results {
    let xnp: f64 = 3.0;
    let tau: f64 = 1.0;
    let xnt: f64 = 193.2;
    let w: f64 = 3.0;

    let mut array_t = Vec::new();
    let mut array_rtmp = Vec::new();

    let mut rt1ic = 500.0;
    while rt1ic <= 40000.0 {
        let vm: f64 = 3000.0;
        let vt: f64 = 1000.0;
        let mut rm1: f64 = 0.0;
        let mut rm2: f64 = 0.0;
        let mut rt1: f64 = rt1ic;
        let mut rt2: f64 = 0.0;
        let mut beta: f64 = 0.0;
        let mut vt1 = -vt * beta.cos();
        let mut vt2 = vt * beta.sin();
        let mut t: f64 = 0.0;

        let mut rtm1 = rt1 - rm1;
        let mut rtm2 = rt2 - rm2;
        let mut rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
        let mut _xlam = rtm2.atan2(rtm1);
        let mut vm1 = vm;
        let mut vm2: f64 = 0.0;
        let vtm1 = vt1 - vm1;
        let vtm2 = vt2 - vm2;
        let mut vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;
        let mut xlamh: f64 = 0.0;
        let mut h: f64 = 0.01;

        while vc > 0.0 {
            if rtm < 1000.0 {
                h = 0.0005;
            }

            let betaold = beta;
            let rt1old = rt1;
            let rt2old = rt2;
            let rm1old = rm1;
            let rm2old = rm2;
            let vm1old = vm1;
            let vm2old = vm2;
            let xlamhold = xlamh;

            // First derivative evaluation
            vt1 = -vt * beta.cos();
            vt2 = vt * beta.sin();
            let betad = xnt * (w * t).sin() / vt;
            rtm1 = rt1 - rm1;
            rtm2 = rt2 - rm2;
            rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
            let vtm1 = vt1 - vm1;
            let vtm2 = vt2 - vm2;
            vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;
            _xlam = rtm2.atan2(rtm1);
            let xlamhd = (_xlam - xlamh) / tau;
            let xnc = xnp * vc * xlamhd;
            let am1 = -xnc * _xlam.sin();
            let am2 = xnc * _xlam.cos();

            // Euler step
            beta += h * betad;
            rt1 += h * vt1;
            rt2 += h * vt2;
            rm1 += h * vm1;
            rm2 += h * vm2;
            vm1 += h * am1;
            vm2 += h * am2;
            xlamh += h * xlamhd;
            t += h;

            // Second derivative for RK2
            vt1 = -vt * beta.cos();
            vt2 = vt * beta.sin();
            let betad = xnt * (w * t).sin() / vt;
            rtm1 = rt1 - rm1;
            rtm2 = rt2 - rm2;
            rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
            let vtm1 = vt1 - vm1;
            let vtm2 = vt2 - vm2;
            vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;
            _xlam = rtm2.atan2(rtm1);
            let xlamhd = (_xlam - xlamh) / tau;
            let xnc = xnp * vc * xlamhd;
            let am1 = -xnc * _xlam.sin();
            let am2 = xnc * _xlam.cos();

            // RK2 averaging
            beta = 0.5 * (betaold + beta + h * betad);
            rt1 = 0.5 * (rt1old + rt1 + h * vt1);
            rt2 = 0.5 * (rt2old + rt2 + h * vt2);
            rm1 = 0.5 * (rm1old + rm1 + h * vm1);
            rm2 = 0.5 * (rm2old + rm2 + h * vm2);
            vm1 = 0.5 * (vm1old + vm1 + h * am1);
            vm2 = 0.5 * (vm2old + vm2 + h * am2);
            xlamh = 0.5 * (xlamhold + xlamh + h * xlamhd);
        }

        let rtmp = if rtm2 > 0.0 { rtm } else { -rtm };

        array_t.push(t);
        array_rtmp.push(rtmp);

        rt1ic += 500.0;
    }

    Results {
        time: array_t,
        miss: array_rtmp,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c29l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.miss.clone(),
    ])?;

    let plot_file = format!("{}/c29l1_miss.png", output_dir);
    let config = PlotConfig::new("Miss Distance vs Flight Time")
        .with_labels("Flight Time (Sec)", "Miss (Ft)");

    let series = vec![
        Series::new(results.time.clone(), results.miss.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C29L1: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c29l1_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
