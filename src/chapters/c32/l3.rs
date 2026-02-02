//! Chapter 32, Lesson 3: Rolling Airframe Guidance (Open Loop)
//!
//! Simulates a rolling airframe missile with constant lateral acceleration
//! using open-loop predicted impact point guidance.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub xm: Vec<f64>,
    pub ym: Vec<f64>,
    pub xt: Vec<f64>,
    pub yt: Vec<f64>,
    pub phideg: Vec<f64>,
}

/// Run the C32L3 simulation
pub fn run() -> Results {
    let xt: f64 = -4000.0;
    let yt: f64 = 5000.0;
    let tf: f64 = 100.0;
    let h: f64 = 0.01;
    let ts: f64 = 0.1;
    let xnc: f64 = 12.0;

    let mut xm: f64 = 0.0;
    let mut ym: f64 = 0.0;
    let mut xmd: f64 = 0.0;
    let mut ymd: f64 = 0.0;
    let mut s: f64 = 0.0;
    let mut t: f64 = 0.0;
    let mut phi: f64 = 0.0;

    let mut array_t = Vec::new();
    let mut array_xm = Vec::new();
    let mut array_ym = Vec::new();
    let mut array_xt = Vec::new();
    let mut array_yt = Vec::new();
    let mut array_phideg = Vec::new();

    while t <= tf {
        let xmold = xm;
        let ymold = ym;
        let xmdold = xmd;
        let ymdold = ymd;

        // First derivative evaluation
        let xmdd = xnc * phi.cos();
        let ymdd = xnc * phi.sin();

        // Euler step
        xm += h * xmd;
        ym += h * ymd;
        xmd += h * xmdd;
        ymd += h * ymdd;
        t += h;

        // Second derivative evaluation
        let xmdd = xnc * phi.cos();
        let ymdd = xnc * phi.sin();

        // RK2 averaging
        xm = (xmold + xm) / 2.0 + 0.5 * h * xmd;
        ym = (ymold + ym) / 2.0 + 0.5 * h * ymd;
        xmd = (xmdold + xmd) / 2.0 + 0.5 * h * xmdd;
        ymd = (ymdold + ymd) / 2.0 + 0.5 * h * ymdd;

        s += h;
        if s >= ts - 0.0001 {
            s = 0.0;
            let tgo = tf - t + 0.0001;
            let rtm1 = xt - xm;
            let rtm2 = yt - ym;
            let _rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
            let vtm1 = -xmd;
            let vtm2 = -ymd;
            phi = (rtm2 + vtm2 * tgo).atan2(rtm1 + vtm1 * tgo);
            let phideg = phi * 57.3;

            array_t.push(t);
            array_xm.push(xm);
            array_ym.push(ym);
            array_xt.push(xt);
            array_yt.push(yt);
            array_phideg.push(phideg);
        }
    }

    Results {
        time: array_t,
        xm: array_xm,
        ym: array_ym,
        xt: array_xt,
        yt: array_yt,
        phideg: array_phideg,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c32l3_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.xm.clone(),
        results.ym.clone(),
        results.xt.clone(),
        results.yt.clone(),
        results.phideg.clone(),
    ])?;

    let plot_file = format!("{}/c32l3_traj.png", output_dir);
    let config = PlotConfig::new("Rolling Airframe Trajectory (Open Loop)")
        .with_labels("X (Ft)", "Y (Ft)");

    let series = vec![
        Series::new(results.xm.clone(), results.ym.clone())
            .with_color(plotters::prelude::BLUE)
            .with_label("Missile"),
        Series::new(results.xt.clone(), results.yt.clone())
            .with_color(plotters::prelude::RED)
            .with_label("Target"),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C32L3: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c32l3_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
