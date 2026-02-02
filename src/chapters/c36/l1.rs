//! Chapter 36, Lesson 1: Impact Angle Control
//!
//! Guidance law with terminal impact angle constraint.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub y: Vec<f64>,
    pub xncg: Vec<f64>,
    pub xlamdeg: Vec<f64>,
}

/// Run the C36L1 simulation
pub fn run() -> Results {
    let xnt: f64 = 0.0;
    let hedeg: f64 = -20.0;
    let xnclim: f64 = 999999.0;
    let pn: i32 = 0;
    let xlamfdeg: f64 = -30.0;
    let vc: f64 = 4000.0;
    let vm: f64 = 3000.0;
    let tf: f64 = 10.0;
    let xnp: f64 = 3.0;

    let xlamf = xlamfdeg / 57.3;
    let mut y: f64 = 0.0;
    let mut yd = -vm * hedeg / 57.3;
    let mut t: f64 = 0.0;
    let h: f64 = 0.001;
    let mut s: f64 = 0.0;

    let mut array_t = Vec::new();
    let mut array_y = Vec::new();
    let mut array_xncg = Vec::new();
    let mut array_xlamdeg = Vec::new();

    while t <= tf - 0.0001 {
        let yold = y;
        let ydold = yd;

        // First derivative evaluation
        let tgo = tf - t + 0.00001;
        let xlam = y / (vc * tgo);
        let xlamd = (y + yd * tgo) / (vc * tgo * tgo);
        let mut xnc = if pn == 1 {
            xnp * vc * xlamd
        } else {
            4.0 * vc * xlamd + xnt + 2.0 * vc * (xlam - xlamf) / tgo
        };
        if xnc > xnclim {
            xnc = xnclim;
        }
        if xnc < -xnclim {
            xnc = -xnclim;
        }
        let ydd = xnt - xnc;

        // Euler step
        y += h * yd;
        yd += h * ydd;
        t += h;

        // Second derivative for RK2
        let tgo = tf - t + 0.00001;
        let xlam = y / (vc * tgo);
        let xlamd = (y + yd * tgo) / (vc * tgo * tgo);
        let mut xnc = if pn == 1 {
            xnp * vc * xlamd
        } else {
            4.0 * vc * xlamd + xnt + 2.0 * vc * (xlam - xlamf) / tgo
        };
        if xnc > xnclim {
            xnc = xnclim;
        }
        if xnc < -xnclim {
            xnc = -xnclim;
        }
        let ydd = xnt - xnc;

        // RK2 averaging
        y = 0.5 * (yold + y + h * yd);
        yd = 0.5 * (ydold + yd + h * ydd);

        s += h;
        if s >= 0.09999 {
            s = 0.0;
            let xlamdeg = xlam * 57.3;
            let xncg = xnc / 32.2;
            array_t.push(t);
            array_y.push(y);
            array_xncg.push(xncg);
            array_xlamdeg.push(xlamdeg);
        }
    }

    Results {
        time: array_t,
        y: array_y,
        xncg: array_xncg,
        xlamdeg: array_xlamdeg,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c36l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.y.clone(),
        results.xncg.clone(),
        results.xlamdeg.clone(),
    ])?;

    let plot_file = format!("{}/c36l1_y.png", output_dir);
    let config = PlotConfig::new("Impact Angle Control - Trajectory")
        .with_labels("Time (Sec)", "Y (Ft)");

    let series = vec![
        Series::new(results.time.clone(), results.y.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C36L1: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c36l1_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
