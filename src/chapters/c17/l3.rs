//! Chapter 17, Lesson 3: Two-Stage Rocket with Lambert Guidance
//!
//! Two-stage rocket trajectory with Lambert guidance to target.
//!
//! MATLAB BUG: Calls `lambertpz` but the function file is `lamberpz.m` and the
//! function is defined as `lambert` not `lambertpz`.
//! Output: 3 columns [T, DISTNM, ALTNM]

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;
use crate::utils::lambert2d::lambert2d;

pub struct Results {
    pub time: Vec<f64>,
    pub distnm: Vec<f64>,
    pub altnm: Vec<f64>,
}

/// Run the C17L3 simulation
pub fn run() -> Results {
    let _left: i32 = 1;
    let mut qboost: bool = true;
    let xisp1: f64 = 300.0;
    let xisp2: f64 = 300.0;
    let xmf1: f64 = 0.90;
    let xmf2: f64 = 0.90;
    let wpay: f64 = 100.0;
    let delv: f64 = 20000.0;
    let delv1 = 0.3333 * delv;
    let delv2 = 0.6667 * delv;
    let amax1: f64 = 20.0;
    let amax2: f64 = 20.0;

    // Stage 2 calculations
    let top2 = wpay * ((delv2 / (xisp2 * 32.2)).exp() - 1.0);
    let bot2 = 1.0 / xmf2 - ((1.0 - xmf2) / xmf2) * (delv2 / (xisp2 * 32.2)).exp();
    let wp2 = top2 / bot2;
    let ws2 = wp2 * (1.0 - xmf2) / xmf2;
    let wtot2 = wp2 + ws2 + wpay;
    let trst2 = amax2 * (wpay + ws2);
    let tb2 = xisp2 * wp2 / trst2;

    // Stage 1 calculations
    let top1 = wtot2 * ((delv1 / (xisp1 * 32.2)).exp() - 1.0);
    let bot1 = 1.0 / xmf1 - ((1.0 - xmf1) / xmf1) * (delv1 / (xisp1 * 32.2)).exp();
    let wp1 = top1 / bot1;
    let ws1 = wp1 * (1.0 - xmf1) / xmf1;
    let wtot = wp1 + ws1 + wtot2;
    let trst1 = amax1 * (wtot2 + ws1);
    let tb1 = xisp1 * wp1 / trst1;

    let h: f64 = 0.01;
    let a: f64 = 2.0926e7;
    let gm: f64 = 1.4077e16;
    let altnm_init: f64 = 0.0;
    let alt_init = altnm_init * 6076.0;
    let angdeg: f64 = 30.0;
    let ang = angdeg / 57.3;
    let _xlongm = ang;

    let mut x = (a + alt_init) * ang.cos();
    let mut y = (a + alt_init) * ang.sin();
    let mut alt = (x * x + y * y).sqrt() - a;
    let xfirst = x;
    let yfirst = y;
    let mut x1: f64 = 0.0;
    let mut y1: f64 = 0.0;
    let mut axt: f64 = 0.0;
    let mut ayt: f64 = 0.0;

    let xlongtdeg: f64 = 45.0;
    let xlongt = xlongtdeg / 57.3;
    let xf = a * xlongt.cos();
    let yf = a * xlongt.sin();
    let tf_target: f64 = 500.0;

    let mut t: f64 = 0.0;
    let mut s: f64 = 0.0;

    let mut array_t = Vec::new();
    let mut array_distnm = Vec::new();
    let mut array_altnm = Vec::new();

    while !(alt < 0.0 && t > 10.0) {
        let x_old = x;
        let y_old = y;
        let x1_old = x1;
        let y1_old = y1;

        // Compute thrust and weight
        let (wgt, trst) = if t < tb1 {
            let w = -wp1 * t / tb1 + wtot;
            (w, trst1)
        } else if t < (tb1 + tb2) {
            let w = -wp2 * t / tb2 + wtot2 + wp2 * tb1 / tb2;
            (w, trst2)
        } else {
            (wpay, 0.0)
        };

        let at = 32.2 * trst / wgt;
        let xd = x1;
        let yd = y1;
        let tembot = (x * x + y * y).powf(1.5);
        let x1d = -gm * x / tembot + axt;
        let y1d = -gm * y / tembot + ayt;

        // Euler step
        x += h * xd;
        y += h * yd;
        x1 += h * x1d;
        y1 += h * y1d;
        t += h;

        // Second derivative
        let xd = x1;
        let yd = y1;
        let tembot = (x * x + y * y).powf(1.5);
        let x1d = -gm * x / tembot + axt;
        let y1d = -gm * y / tembot + ayt;

        // RK2 averaging
        x = (x_old + x) / 2.0 + 0.5 * h * xd;
        y = (y_old + y) / 2.0 + 0.5 * h * yd;
        x1 = (x1_old + x1) / 2.0 + 0.5 * h * x1d;
        y1 = (y1_old + y1) / 2.0 + 0.5 * h * y1d;
        alt = (x * x + y * y).sqrt() - a;

        if qboost {
            let tgolam = tf_target - t;
            let xlongm = y.atan2(x);
            let result = lambert2d(x, y, tgolam, xf, yf, xlongm, xlongt);
            let vrx = result.vrx;
            let vry = result.vry;
            let delx = vrx - x1;
            let dely = vry - y1;
            let del = (delx * delx + dely * dely).sqrt();

            if trst > 0.0 && del > 500.0 {
                axt = at * delx / del;
                ayt = at * dely / del;
            } else if del < 500.0 {
                qboost = false;
                axt = 0.0;
                ayt = 0.0;
                x1 = vrx;
                y1 = vry;
            } else {
                qboost = false;
                axt = 0.0;
                ayt = 0.0;
            }
        }

        s += h;
        if s >= 9.99999 {
            s = 0.0;
            let r = (x * x + y * y).sqrt();
            let rf = (xfirst * xfirst + yfirst * yfirst).sqrt();
            let cbeta = (x * xfirst + y * yfirst) / (r * rf);
            let beta = cbeta.acos();
            let distnm = a * beta / 6076.0;
            let altnm = ((x * x + y * y).sqrt() - a) / 6076.0;

            array_t.push(t);
            array_distnm.push(distnm);
            array_altnm.push(altnm);
        }
    }

    Results {
        time: array_t,
        distnm: array_distnm,
        altnm: array_altnm,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c17l3_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.distnm.clone(),
        results.altnm.clone(),
    ])?;

    let plot_file = format!("{}/c17l3_trajectory.png", output_dir);
    let config = PlotConfig::new("Lambert Guided Trajectory")
        .with_labels("Downrange (Nmi)", "Altitude (Nmi)");
    let series = vec![
        Series::new(results.distnm.clone(), results.altnm.clone())
            .with_color(plotters::prelude::BLUE),
    ];
    line_plot(&plot_file, &config, &series).ok();

    println!("C17L3: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c17l3_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
