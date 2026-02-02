//! Chapter 15, Lesson 2: Flat Earth vs Spherical Earth
//!
//! Comparison of polar and Cartesian equations for spherical Earth.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub spolarnm: Vec<f64>,
    pub altpolarnm: Vec<f64>,
    pub distnm: Vec<f64>,
    pub altnm: Vec<f64>,
}

/// Run the C15L2 simulation
pub fn run() -> Results {
    let h: f64 = 0.01;
    let a: f64 = 2.0926e7;  // Earth radius in ft
    let gm: f64 = 1.4077e16; // Gravitational parameter
    let gam: f64 = 45.0;    // Initial flight path angle (deg)
    let altnm_init: f64 = 0.0;
    let v: f64 = 24000.0;   // Initial velocity (ft/s)
    let angdeg: f64 = 0.0;
    let ang = angdeg / 57.3;

    let vrx = v * (1.5708 - gam / 57.3 + ang).cos();
    let vry = v * (1.5708 - gam / 57.3 + ang).sin();
    let alt = altnm_init / 6076.0;
    let mut s: f64 = 0.0;

    // Polar coordinates
    let mut r0 = a + alt;
    let mut r1 = v * (gam / 57.3).sin();
    let mut psi: f64 = 0.0;

    // Cartesian coordinates
    let mut x = (a + alt) * ang.cos();
    let mut y = (a + alt) * ang.sin();
    let xfirst = x;
    let yfirst = y;
    let mut x1 = vrx;
    let mut y1 = vry;

    let mut t: f64 = 0.0;
    let mut altnm = altnm_init;

    let mut array_t = Vec::new();
    let mut array_spolarnm = Vec::new();
    let mut array_altpolarnm = Vec::new();
    let mut array_distnm = Vec::new();
    let mut array_altnm = Vec::new();

    while altnm > -0.0001 {
        let r0_old = r0;
        let r1_old = r1;
        let psi_old = psi;
        let x_old = x;
        let y_old = y;
        let x1_old = x1;
        let y1_old = y1;

        // First derivative evaluation
        let psid = (a + alt) * v * (gam / 57.3).cos() / (r0 * r0);
        let r1d = -gm / (r0 * r0) + r0 * psid * psid;
        let r0d = r1;
        let tembot = (x * x + y * y).powf(1.5);
        let x1d = -gm * x / tembot;
        let y1d = -gm * y / tembot;
        let xd = x1;
        let yd = y1;

        // Euler step
        r0 += h * r0d;
        r1 += h * r1d;
        psi += h * psid;
        x += h * xd;
        y += h * yd;
        x1 += h * x1d;
        y1 += h * y1d;
        t += h;

        // Second derivative
        let psid = (a + alt) * v * (gam / 57.3).cos() / (r0 * r0);
        let r1d = -gm / (r0 * r0) + r0 * psid * psid;
        let r0d = r1;
        let tembot = (x * x + y * y).powf(1.5);
        let x1d = -gm * x / tembot;
        let y1d = -gm * y / tembot;
        let xd = x1;
        let yd = y1;

        // RK2 averaging
        r0 = (r0_old + r0) / 2.0 + 0.5 * h * r0d;
        r1 = (r1_old + r1) / 2.0 + 0.5 * h * r1d;
        psi = (psi_old + psi) / 2.0 + 0.5 * h * psid;
        x = (x_old + x) / 2.0 + 0.5 * h * xd;
        y = (y_old + y) / 2.0 + 0.5 * h * yd;
        x1 = (x1_old + x1) / 2.0 + 0.5 * h * x1d;
        y1 = (y1_old + y1) / 2.0 + 0.5 * h * y1d;

        s += h;
        if s >= 9.99999 {
            s = 0.0;
            let spolarnm = a * psi / 6076.0;
            let altpolarnm = (r0 - a) / 6076.0;
            altnm = ((x * x + y * y).sqrt() - a) / 6076.0;
            let r = (x * x + y * y).sqrt();
            let rf = (xfirst * xfirst + yfirst * yfirst).sqrt();
            let cbeta = (x * xfirst + y * yfirst) / (r * rf);
            let beta = cbeta.acos();
            let distnm = a * beta / 6076.0;

            array_t.push(t);
            array_spolarnm.push(spolarnm);
            array_altpolarnm.push(altpolarnm);
            array_distnm.push(distnm);
            array_altnm.push(altnm);
        }
    }

    Results {
        time: array_t,
        spolarnm: array_spolarnm,
        altpolarnm: array_altpolarnm,
        distnm: array_distnm,
        altnm: array_altnm,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c15l2_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.spolarnm.clone(),
        results.altpolarnm.clone(),
        results.distnm.clone(),
        results.altnm.clone(),
    ])?;

    // Trajectory comparison plot
    let plot_file = format!("{}/c15l2_trajectory.png", output_dir);
    let config = PlotConfig::new("Polar vs Cartesian Trajectory")
        .with_labels("Downrange (Nmi)", "Altitude (Nmi)");
    let series = vec![
        Series::new(results.spolarnm.clone(), results.altpolarnm.clone())
            .with_label("Polar")
            .with_color(plotters::prelude::BLUE),
        Series::new(results.distnm.clone(), results.altnm.clone())
            .with_label("Cartesian")
            .with_color(plotters::prelude::RED),
    ];
    line_plot(&plot_file, &config, &series).ok();

    println!("C15L2: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c15l2_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
