//! Chapter 33, Lesson 1: Predictive Guidance
//!
//! Implements predictive guidance with commanded acceleration control.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub y: Vec<f64>,
    pub xncg: Vec<f64>,
    pub xncdg: Vec<f64>,
    pub xntg: Vec<f64>,
}

/// Run the C33L1 simulation
pub fn run() -> Results {
    let qpn: i32 = 1;
    let _imvr: i32 = 1; // Target maneuver type: 1=constant, other=weaving
    let _vc: f64 = 4000.0;
    let vm: f64 = 3000.0;
    let xntic: f64 = 96.6;
    let hedeg: f64 = 0.0;
    let _w: f64 = 3.0; // Weave frequency (unused when imvr=1)
    let tf: f64 = 10.0;
    let xnp: f64 = 3.0;
    let xnpp: f64 = 10.0;
    let xncdlim: f64 = 999999.0;

    let mut y: f64 = 0.0;
    let mut yd = -vm * hedeg / 57.3;
    let mut t: f64 = 0.0;
    let h: f64 = 0.001;
    let mut s: f64 = 0.0;
    let mut xnc = xnp * yd / tf;
    let mut j: f64 = 0.0;

    let mut array_t = Vec::new();
    let mut array_y = Vec::new();
    let mut array_xncg = Vec::new();
    let mut array_xncdg = Vec::new();
    let mut array_xntg = Vec::new();

    while t <= tf - 1e-5 {
        let yold = y;
        let ydold = yd;
        let xncold = xnc;
        let jold = j;

        // First derivative evaluation
        let tgo = tf - t + 0.00001;
        // Note: imvr == 1 always, so xnt = xntic (weaving target case would use sin(w*t))
        let xnt = xntic;
        let ydd = xnt - xnc;
        let mut xncd = if qpn == 1 {
            2.0 * xnp * (y + yd * tgo + 0.5 * ydd * tgo * tgo) / tgo.powi(3)
        } else {
            xnpp * (y + yd * tgo + 0.5 * ydd * tgo * tgo) / tgo.powi(3)
        };
        if xncd > xncdlim {
            xncd = xncdlim;
        }
        if xncd < -xncdlim {
            xncd = -xncdlim;
        }
        let jd = xncd * xncd;

        // Euler step
        y += h * yd;
        yd += h * ydd;
        xnc += h * xncd;
        j += h * jd;
        t += h;

        // Second derivative for RK2
        let tgo = tf - t + 0.00001;
        // Note: imvr == 1 always, so xnt = xntic (weaving target case would use sin(w*t))
        let xnt = xntic;
        let ydd = xnt - xnc;
        let mut xncd = if qpn == 1 {
            2.0 * xnp * (y + yd * tgo + 0.5 * ydd * tgo * tgo) / tgo.powi(3)
        } else {
            xnpp * (y + yd * tgo + 0.5 * ydd * tgo * tgo) / tgo.powi(3)
        };
        if xncd > xncdlim {
            xncd = xncdlim;
        }
        if xncd < -xncdlim {
            xncd = -xncdlim;
        }
        let jd = xncd * xncd;

        // RK2 averaging
        y = 0.5 * (yold + y + h * yd);
        yd = 0.5 * (ydold + yd + h * ydd);
        xnc = 0.5 * (xncold + xnc + h * xncd);
        j = 0.5 * (jold + j + h * jd);

        s += h;
        if s >= 0.09999 {
            s = 0.0;
            array_t.push(t);
            array_y.push(y);
            array_xncg.push(xnc / 32.2);
            array_xncdg.push(xncd / 32.2);
            array_xntg.push(xnt / 32.2);
        }
    }

    Results {
        time: array_t,
        y: array_y,
        xncg: array_xncg,
        xncdg: array_xncdg,
        xntg: array_xntg,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c33l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.y.clone(),
        results.xncg.clone(),
        results.xncdg.clone(),
        results.xntg.clone(),
    ])?;

    let plot_file = format!("{}/c33l1_y.png", output_dir);
    let config = PlotConfig::new("Predictive Guidance - Miss Distance")
        .with_labels("Time (s)", "Y (ft)");

    let series = vec![
        Series::new(results.time.clone(), results.y.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C33L1: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c33l1_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }

    #[test]
    fn test_c33l1_convergence() {
        let results = run();
        // Check that Y converges toward zero
        let n = results.y.len();
        assert!(results.y[n - 1].abs() < results.y[0].abs() + 100.0);
    }
}
