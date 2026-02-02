//! Chapter 25, Lesson 3: Flexible Body Effects with Two Modes and Notch Filters
//!
//! Simulates flexible body with two modes and optional notch filters.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub thdtot: Vec<f64>,  // Total body rate
}

/// Run the C25L3 simulation
pub fn run() -> Results {
    let ts: f64 = 0.001;
    let xk3: f64 = -0.362;
    let ta: f64 = 0.831;
    let b11: f64 = 0.00461;
    let b12: f64 = 0.00136;
    let xkfb: f64 = 0.00134;
    let wzfb: f64 = 395.0;
    let z1: f64 = 0.015;
    let w1: f64 = 259.0;
    let xkfb2: f64 = 0.000664;
    let wzfb2: f64 = 255.0;
    let z2: f64 = 0.022;
    let w2: f64 = 649.0;
    let wt: f64 = w1;
    let zt: f64 = z1;
    let wt2: f64 = w2;
    let zt2: f64 = z2;
    let wb: f64 = w1;
    let zb: f64 = 0.5;
    let wb2: f64 = w2;
    let zb2: f64 = 0.5;
    let tf: f64 = 1.0;
    let xin: f64 = 1.0;
    // WACT=100.; then overwritten to 400
    let wact: f64 = 400.0;
    let zact: f64 = 0.7;
    // XKR=.1 then .3
    let xkr: f64 = 0.3;
    // NOTCH=3 (no filters), NOTCH=1 (first mode filter), NOTCH=2 (both filters)
    let notch: i32 = 2;
    let xk: f64 = (1.0 - xkr * xk3) / (xkr * xk3);
    // FB1=1. (First mode) FB2=1 (Second mode)
    let fb1: f64 = 1.0;
    let fb2: f64 = 1.0;

    let mut t: f64 = 0.0;
    let mut s: f64 = 0.0;
    let h: f64 = 0.00001;
    let mut e: f64 = 0.0;
    let mut ed: f64 = 0.0;
    let mut e1: f64 = 0.0;
    let mut e1d: f64 = 0.0;
    let mut e1dd: f64 = 0.0;
    let mut del: f64 = 0.0;
    let mut deld: f64 = 0.0;
    let mut delc: f64 = 0.0;
    let mut e2: f64 = 0.0;
    let mut e2d: f64 = 0.0;
    let mut e3: f64 = 0.0;
    let mut e3d: f64 = 0.0;
    let mut e3dd: f64 = 0.0;
    let mut e4: f64 = 0.0;
    let mut e4d: f64 = 0.0;
    let mut thdfil1: f64 = 0.0;
    let mut thdfil: f64 = 0.0;

    let mut array_t = Vec::new();
    let mut array_thdtot = Vec::new();

    while t <= tf {
        let eold = e;
        let edold = ed;
        let e1old = e1;
        let e1dold = e1d;
        let e1ddold = e1dd;
        let delold = del;
        let deldold = deld;
        let e2old = e2;
        let e2dold = e2d;
        let e3old = e3;
        let e3dold = e3d;
        let e3ddold = e3dd;
        let e4old = e4;
        let e4dold = e4d;

        // First derivative evaluation
        let deldd = wact * wact * (delc - del - 2.0 * zact * deld / wact);
        let edd = (del - e - b11 * ed) / b12;
        let e1ddd = w1 * w1 * (deld - e1d - 2.0 * z1 * e1dd / w1);
        let e3ddd = w2 * w2 * (deld - e3d - 2.0 * z2 * e3dd / w2);
        let thdrb = xk3 * (e + ta * ed);
        let thdfb = xkfb * (e1d + e1ddd / (wzfb * wzfb));
        let thdfb2 = xkfb2 * (e3d + e3ddd / (wzfb2 * wzfb2));
        let thdtot = thdrb + fb1 * thdfb + fb2 * thdfb2;

        if notch == 1 {
            delc = xkr * (xk * xin + thdfil1);
        } else if notch == 2 {
            delc = xkr * (xk * xin + thdfil);
        } else {
            delc = xkr * (xk * xin + thdtot);
        }

        let e2dd = wb * wb * (thdtot - e2 - 2.0 * zb * e2d / wb);
        thdfil1 = e2 + 2.0 * zt * e2d / wt + e2dd / (wt * wt);
        let e4dd = wb2 * wb2 * (thdfil1 - e4 - 2.0 * zb2 * e4d / wb2);
        thdfil = e4 + 2.0 * zt2 * e4d / wt2 + e4dd / (wt2 * wt2);

        // Euler step
        e += h * ed;
        ed += h * edd;
        e1 += h * e1d;
        e1d += h * e1dd;
        e1dd += h * e1ddd;
        del += h * deld;
        deld += h * deldd;
        e2 += h * e2d;
        e2d += h * e2dd;
        e3 += h * e3d;
        e3d += h * e3dd;
        e3dd += h * e3ddd;
        e4 += h * e4d;
        e4d += h * e4dd;
        t += h;

        // Second derivative for RK2
        let deldd = wact * wact * (delc - del - 2.0 * zact * deld / wact);
        let edd = (del - e - b11 * ed) / b12;
        let e1ddd = w1 * w1 * (deld - e1d - 2.0 * z1 * e1dd / w1);
        let e3ddd = w2 * w2 * (deld - e3d - 2.0 * z2 * e3dd / w2);
        let thdrb = xk3 * (e + ta * ed);
        let thdfb = xkfb * (e1d + e1ddd / (wzfb * wzfb));
        let thdfb2 = xkfb2 * (e3d + e3ddd / (wzfb2 * wzfb2));
        let thdtot = thdrb + fb1 * thdfb + fb2 * thdfb2;

        if notch == 1 {
            delc = xkr * (xk * xin + thdfil1);
        } else if notch == 2 {
            delc = xkr * (xk * xin + thdfil);
        } else {
            delc = xkr * (xk * xin + thdtot);
        }

        let e2dd = wb * wb * (thdtot - e2 - 2.0 * zb * e2d / wb);
        thdfil1 = e2 + 2.0 * zt * e2d / wt + e2dd / (wt * wt);
        let e4dd = wb2 * wb2 * (thdfil1 - e4 - 2.0 * zb2 * e4d / wb2);
        thdfil = e4 + 2.0 * zt2 * e4d / wt2 + e4dd / (wt2 * wt2);

        // RK2 averaging
        e = 0.5 * (eold + e + h * ed);
        ed = 0.5 * (edold + ed + h * edd);
        e1 = 0.5 * (e1old + e1 + h * e1d);
        e1d = 0.5 * (e1dold + e1d + h * e1dd);
        e1dd = 0.5 * (e1ddold + e1dd + h * e1ddd);
        del = 0.5 * (delold + del + h * deld);
        deld = 0.5 * (deldold + deld + h * deldd);
        e2 = 0.5 * (e2old + e2 + h * e2d);
        e2d = 0.5 * (e2dold + e2d + h * e2dd);
        e3 = 0.5 * (e3old + e3 + h * e3d);
        e3d = 0.5 * (e3dold + e3d + h * e3dd);
        e3dd = 0.5 * (e3ddold + e3dd + h * e3ddd);
        e4 = 0.5 * (e4old + e4 + h * e4d);
        e4d = 0.5 * (e4dold + e4d + h * e4dd);

        s += h;
        if s >= ts - 0.00001 {
            s = 0.0;
            array_t.push(t);
            array_thdtot.push(thdtot);
        }
    }

    Results {
        time: array_t,
        thdtot: array_thdtot,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c25l3_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.thdtot.clone(),
    ])?;

    let plot_file = format!("{}/c25l3_thdtot.png", output_dir);
    let config = PlotConfig::new("Flexible Body - Two Modes with Notch Filters")
        .with_labels("Time (Sec)", "Total Body Rate (deg/s)");

    let series = vec![
        Series::new(results.time.clone(), results.thdtot.clone())
            .with_color(plotters::prelude::BLUE)
            .with_label("Total Body Rate"),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C25L3: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c25l3_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
