//! Chapter 17, Lesson 6: Ballistic Trajectory Propagation
//!
//! Propagate ballistic trajectory with given initial conditions.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub distnm: Vec<f64>,
    pub altnm: Vec<f64>,
}

/// Run the C17L6 simulation
pub fn run() -> Results {
    let pi: f64 = std::f64::consts::PI;
    let degrad = 360.0 / (2.0 * pi);
    let h: f64 = 0.01;
    let a: f64 = 2.0926e7;
    let gm: f64 = 1.4077e16;
    let gam: f64 = 30.0;
    let altnm_init: f64 = 0.0;
    let v: f64 = 24000.0;
    let tf: f64 = 2200.0;
    let alt = altnm_init * 6076.0;
    let ang: f64 = 0.0;

    let vrx = v * (pi / 2.0 - gam / degrad + ang).cos();
    let vry = v * (pi / 2.0 - gam / degrad + ang).sin();

    let mut x = (a + alt) * ang.cos();
    let mut y = (a + alt) * ang.sin();
    let xfirst = x;
    let yfirst = y;
    let mut x1 = vrx;
    let mut y1 = vry;
    let mut t: f64 = 0.0;
    let mut scount: f64 = 0.0;

    let mut array_t = Vec::new();
    let mut array_distnm = Vec::new();
    let mut array_altnm = Vec::new();

    while t < (tf - 0.0001) {
        let x_old = x;
        let y_old = y;
        let x1_old = x1;
        let y1_old = y1;

        // First derivative
        let tembot = (x * x + y * y).powf(1.5);
        let x1d = -gm * x / tembot;
        let y1d = -gm * y / tembot;
        let xd = x1;
        let yd = y1;

        // Euler step
        x += h * xd;
        y += h * yd;
        x1 += h * x1d;
        y1 += h * y1d;
        t += h;

        // Second derivative
        let tembot = (x * x + y * y).powf(1.5);
        let x1d = -gm * x / tembot;
        let y1d = -gm * y / tembot;
        let xd = x1;
        let yd = y1;

        // RK2 averaging
        x = (x_old + x) / 2.0 + 0.5 * h * xd;
        y = (y_old + y) / 2.0 + 0.5 * h * yd;
        x1 = (x1_old + x1) / 2.0 + 0.5 * h * x1d;
        y1 = (y1_old + y1) / 2.0 + 0.5 * h * y1d;

        scount += h;
        if scount >= 1.99999 {
            scount = 0.0;
            let altnm = ((x * x + y * y).sqrt() - a) / 6076.0;
            let r = (x * x + y * y).sqrt();
            let rf = (xfirst * xfirst + yfirst * yfirst).sqrt();
            let cbeta = (x * xfirst + y * yfirst) / (r * rf);
            let beta = cbeta.acos();
            let distnm = a * beta / 6076.0;

            array_t.push(t);
            array_distnm.push(distnm);
            array_altnm.push(altnm);
        }
    }

    println!("X = {:.6e}", x);
    println!("Y = {:.6e}", y);

    Results {
        time: array_t,
        distnm: array_distnm,
        altnm: array_altnm,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c17l6_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.distnm.clone(),
        results.altnm.clone(),
    ])?;

    let plot_file = format!("{}/c17l6_trajectory.png", output_dir);
    let config = PlotConfig::new("Ballistic Trajectory")
        .with_labels("Downrange (Nmi)", "Altitude (Nmi)");
    let series = vec![
        Series::new(results.distnm.clone(), results.altnm.clone())
            .with_color(plotters::prelude::BLUE),
    ];
    line_plot(&plot_file, &config, &series).ok();

    println!("C17L6: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c17l6_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
