//! Chapter 45, Lesson 1: KKV vs Accelerating Target
//!
//! Simulates a Kill Vehicle (KKV) intercepting an accelerating target.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub xnt: Vec<f64>,
    pub delv: Vec<f64>,
}

/// Run the C45L1 simulation
pub fn run() -> Results {
    let ioption: i32 = 0;
    let itgt: i32 = 1;  // Target type: 1 = theater missile
    let xntav: f64 = 117.6;

    let tf = if itgt == 1 { 180.0 } else { 240.0 };

    let pred: f64 = 0.0 * 3280.0;
    let vm: f64 = 9000.0;
    let vc_init: f64 = 18000.0;
    let xntmax: f64 = 9.0 * 32.2;
    let xncmax: f64 = 966.0;
    let apn: f64 = 0.0;
    let hedeg = -57.3 * pred / (vm * tf);
    let xnp: f64 = 3.0;
    let h: f64 = 0.001;

    let mut yd = -vm * hedeg / 57.3;
    let mut y: f64 = 0.0;
    let mut t: f64 = 0.0;
    let mut s: f64 = 0.0;
    let mut delv: f64 = 0.0;
    let mut _sum: f64 = 0.0;
    let mut _xn: f64 = 0.0;

    let mut array_t = Vec::new();
    let mut array_xnt = Vec::new();
    let mut array_delv = Vec::new();

    while t < tf - 0.0001 {
        let yold = y;
        let ydold = yd;
        let delvold = delv;

        // First derivative evaluation
        let tgo = tf - t + 0.00001;

        // Target acceleration
        let xnt = if itgt == 1 {
            if ioption == 0 {
                xntmax * (t / tf).powi(2)
            } else if t < 180.0 {
                let wgt = -212.0 * t + 44000.0;
                let trst = 54100.0;
                32.2 * trst / wgt
            } else {
                let wgt = 3300.0;
                let trst = 0.0;
                32.2 * trst / wgt
            }
        } else if ioption == 0 {
            xntav
        } else if t < 120.0 {
            let wgt = -2622.0 * t + 440660.0;
            let trst = 725850.0;
            32.2 * trst / wgt
        } else if t < 240.0 {
            let wgt = -642.0 * t + 168120.0;
            let trst = 182250.0;
            32.2 * trst / wgt
        } else {
            let wgt = 5500.0;
            let trst = 0.0;
            32.2 * trst / wgt
        };

        let xlamd = (y + yd * tgo) / (vc_init * tgo * tgo);
        let mut xnc = xnp * vc_init * xlamd + 0.5 * apn * xnp * xnt;

        if xnc > xncmax {
            xnc = xncmax;
        }
        if xnc < -xncmax {
            xnc = -xncmax;
        }

        let delvd = xnc.abs();
        let ydd = xnt - xnc;

        // Euler step
        y += h * yd;
        yd += h * ydd;
        delv += h * delvd;
        t += h;

        // Second derivative evaluation
        let tgo2 = tf - t + 0.00001;

        let xnt2 = if itgt == 1 {
            if ioption == 0 {
                xntmax * (t / tf).powi(2)
            } else if t < 180.0 {
                let wgt = -212.0 * t + 44000.0;
                let trst = 54100.0;
                32.2 * trst / wgt
            } else {
                let wgt = 3300.0;
                let trst = 0.0;
                32.2 * trst / wgt
            }
        } else if ioption == 0 {
            xntav
        } else if t < 120.0 {
            let wgt = -2622.0 * t + 440660.0;
            let trst = 725850.0;
            32.2 * trst / wgt
        } else if t < 240.0 {
            let wgt = -642.0 * t + 168120.0;
            let trst = 182250.0;
            32.2 * trst / wgt
        } else {
            let wgt = 5500.0;
            let trst = 0.0;
            32.2 * trst / wgt
        };

        let xlamd2 = (y + yd * tgo2) / (vc_init * tgo2 * tgo2);
        let mut xnc2 = xnp * vc_init * xlamd2 + 0.5 * apn * xnp * xnt2;

        if xnc2 > xncmax {
            xnc2 = xncmax;
        }
        if xnc2 < -xncmax {
            xnc2 = -xncmax;
        }

        let delvd2 = xnc2.abs();
        let ydd2 = xnt2 - xnc2;

        // RK2 averaging
        y = 0.5 * (yold + y + h * yd);
        yd = 0.5 * (ydold + yd + h * ydd2);
        delv = 0.5 * (delvold + delv + h * delvd2);

        s += h;

        if s >= 0.09999 {
            s = 0.0;
            _sum += xnt;
            _xn += 1.0;

            array_t.push(t);
            array_xnt.push(xnt / 32.2);
            array_delv.push(delv / 3.28);
        }
    }

    Results {
        time: array_t,
        xnt: array_xnt,
        delv: array_delv,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c45l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.xnt.clone(),
        results.delv.clone(),
    ])?;

    let plot_file = format!("{}/c45l1_acceleration.png", output_dir);
    let config = PlotConfig::new("KKV vs Accelerating Target")
        .with_labels("Missile Flight Time (s)", "Acceleration (g)");

    let series = vec![
        Series::new(results.time.clone(), results.xnt.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C45L1: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c45l1_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
