//! Chapter 15, Lesson 1: Flat Earth vs Spherical Earth
//!
//! Compares flat earth and spherical earth trajectory models
//! for ballistic flight.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub rt1nm: Vec<f64>,   // Flat earth downrange (Nmi)
    pub rt2nm: Vec<f64>,   // Flat earth altitude (Nmi)
    pub distnm: Vec<f64>,  // Spherical earth great circle distance (Nmi)
    pub altnm: Vec<f64>,   // Spherical earth altitude (Nmi)
}

/// Run the C15L1 simulation
pub fn run() -> Results {
    let h: f64 = 0.01;
    let a: f64 = 2.0926e7;    // Earth radius (ft)
    let gm: f64 = 1.4077e16;  // Gravitational parameter
    let gam: f64 = 45.0;      // Launch angle (deg)
    let v: f64 = 3000.0;      // Initial velocity (ft/s)
    let g: f64 = 32.2;

    let altnm_init: f64 = 0.0;
    let alt = altnm_init / 6076.0;
    let ang: f64 = 0.0;

    // Initial velocities
    let vrx = v * (1.5708 - gam / 57.3 + ang).cos();
    let vry = v * (1.5708 - gam / 57.3 + ang).sin();

    let mut t: f64 = 0.0;
    let mut s: f64 = 0.0;

    // Flat earth state
    let mut rt1 = alt * ang.cos();
    let mut rt2 = alt * ang.sin();
    let mut vt1 = vrx;
    let mut vt2 = vry;

    // Spherical earth state
    let mut x = (a + alt) * ang.cos();
    let mut y = (a + alt) * ang.sin();
    let xfirst = x;
    let yfirst = y;
    let mut x1 = vrx;
    let mut y1 = vry;

    let mut altnm: f64 = altnm_init;

    let mut array_t = Vec::new();
    let mut array_rt1nm = Vec::new();
    let mut array_rt2nm = Vec::new();
    let mut array_distnm = Vec::new();
    let mut array_altnm = Vec::new();

    while altnm > -0.0001 {
        let rt1old = rt1;
        let rt2old = rt2;
        let vt1old = vt1;
        let vt2old = vt2;
        let xold = x;
        let yold = y;
        let x1old = x1;
        let y1old = y1;

        // First derivative evaluation
        // Flat earth
        let at1: f64 = 0.0;
        let at2 = -g;
        let rt1d = vt1;
        let rt2d = vt2;
        let vt1d = at1;
        let vt2d = at2;

        // Spherical earth
        let tembot = (x * x + y * y).powf(1.5);
        let x1d = -gm * x / tembot;
        let y1d = -gm * y / tembot;
        let xd = x1;
        let yd = y1;

        // Euler step
        rt1 += h * rt1d;
        rt2 += h * rt2d;
        vt1 += h * vt1d;
        vt2 += h * vt2d;
        x += h * xd;
        y += h * yd;
        x1 += h * x1d;
        y1 += h * y1d;
        t += h;

        // Second derivative for RK2
        // Flat earth
        let at1: f64 = 0.0;
        let at2 = -g;
        let rt1d = vt1;
        let rt2d = vt2;
        let vt1d = at1;
        let vt2d = at2;

        // Spherical earth
        let tembot = (x * x + y * y).powf(1.5);
        let x1d = -gm * x / tembot;
        let y1d = -gm * y / tembot;
        let xd = x1;
        let yd = y1;

        // RK2 averaging
        rt1 = (rt1old + rt1) / 2.0 + 0.5 * h * rt1d;
        rt2 = (rt2old + rt2) / 2.0 + 0.5 * h * rt2d;
        vt1 = (vt1old + vt1) / 2.0 + 0.5 * h * vt1d;
        vt2 = (vt2old + vt2) / 2.0 + 0.5 * h * vt2d;
        x = (xold + x) / 2.0 + 0.5 * h * xd;
        y = (yold + y) / 2.0 + 0.5 * h * yd;
        x1 = (x1old + x1) / 2.0 + 0.5 * h * x1d;
        y1 = (y1old + y1) / 2.0 + 0.5 * h * y1d;

        s += h;

        if s >= 1.99999 {
            s = 0.0;
            let rt1nm = rt1 / 6076.0;
            let rt2nm = rt2 / 6076.0;
            altnm = ((x * x + y * y).sqrt() - a) / 6076.0;

            let r = (x * x + y * y).sqrt();
            let rf = (xfirst * xfirst + yfirst * yfirst).sqrt();
            let cbeta = (x * xfirst + y * yfirst) / (r * rf);
            let beta = cbeta.acos();
            let distnm = a * beta / 6076.0;

            array_t.push(t);
            array_rt1nm.push(rt1nm);
            array_rt2nm.push(rt2nm);
            array_distnm.push(distnm);
            array_altnm.push(altnm);
        }
    }

    Results {
        time: array_t,
        rt1nm: array_rt1nm,
        rt2nm: array_rt2nm,
        distnm: array_distnm,
        altnm: array_altnm,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c15l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.rt1nm.clone(),
        results.rt2nm.clone(),
        results.distnm.clone(),
        results.altnm.clone(),
    ])?;

    let plot_file = format!("{}/c15l1_plot.png", output_dir);
    let config = PlotConfig::new("Flat Earth vs Spherical Earth")
        .with_labels("Downrange (Nmi)", "Altitude (Nmi)");

    let series = vec![
        Series::new(results.rt1nm.clone(), results.rt2nm.clone())
            .with_color(plotters::prelude::BLUE)
            .with_label("Flat Earth"),
        Series::new(results.distnm.clone(), results.altnm.clone())
            .with_color(plotters::prelude::RED)
            .with_label("Spherical Earth"),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C15L1: Simulation finished");
    println!("  Data saved to: {}", data_file);
    println!("  Plot saved to: {}", plot_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c15l1_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }

    #[test]
    fn test_c15l1_trajectories_similar() {
        let results = run();
        // For short ranges, flat and spherical should be similar
        if !results.rt1nm.is_empty() && !results.distnm.is_empty() {
            let diff = (results.rt1nm[0] - results.distnm[0]).abs();
            assert!(diff < 1.0);  // Within 1 nautical mile at first sample
        }
    }
}
