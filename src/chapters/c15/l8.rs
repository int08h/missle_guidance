//! Chapter 15, Lesson 8: High Trajectory (Lob) Ballistic Flight
//!
//! Ballistic trajectory with high (170 deg) flight path angle.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub distrtkm: Vec<f64>,
    pub altkm: Vec<f64>,
    pub vkm: Vec<f64>,
}

/// Run the C15L8 simulation
pub fn run() -> Results {
    let h: f64 = 0.01;
    let gm: f64 = 1.4077e16;
    let a: f64 = 2.0926e7;
    let gamdeg: f64 = 170.0;
    let gam = gamdeg / 57.3;
    let distkm: f64 = 10000.0;
    let angdeg: f64 = 0.0;
    let ang = angdeg / 57.3;
    let phi = distkm * 3280.0 / a;
    let altkm_init: f64 = 0.0;
    let mut alt = altkm_init * 3280.0;
    let r0 = a + alt;

    let top = gm * (1.0 - phi.cos());
    let temp = r0 * gam.cos() / a - (phi + gam).cos();
    let bot = r0 * gam.cos() * temp;
    let v = (top / bot).sqrt();

    let vrx = v * (1.5708 - gam + ang).cos();
    let vry = v * (1.5708 - gam + ang).sin();

    let mut x = (a + alt) * ang.cos();
    let mut y = (a + alt) * ang.sin();
    let xfirst = x;
    let yfirst = y;
    let mut x1 = vrx;
    let mut y1 = vry;
    let mut t: f64 = 0.0;
    let mut scount: f64 = 0.0;
    let mut ipz: i32 = 1;
    let mut qfirst: bool = true;
    let mut betaold: f64 = 0.0;

    let mut array_distrtkm = Vec::new();
    let mut array_altkm = Vec::new();
    let mut array_vkm = Vec::new();

    while !(t > 10.0 && alt < 0.0) {
        let _x_old = x;
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
        let _xd = x1;
        let yd = y1;

        // RK2 averaging (note: MATLAB has bug, missing X update)
        y = (y_old + y) / 2.0 + 0.5 * h * yd;
        x1 = (x1_old + x1) / 2.0 + 0.5 * h * x1d;
        y1 = (y1_old + y1) / 2.0 + 0.5 * h * y1d;
        alt = (x * x + y * y).sqrt() - a;

        scount += h;
        if scount >= 9.9999 {
            scount = 0.0;
            let altkm = alt / 3280.0;
            let vkm = (x1 * x1 + y1 * y1).sqrt() / 3280.0;
            let r = (x * x + y * y).sqrt();
            let rf = (xfirst * xfirst + yfirst * yfirst).sqrt();
            let mut cbeta = (x * xfirst + y * yfirst) / (r * rf);
            if cbeta > 1.0 {
                cbeta = 1.0;
            }
            let mut beta = cbeta.acos();
            if ipz == 0 {
                beta = 2.0 * std::f64::consts::PI - beta;
            }
            let distrtkm = a * beta / 3280.0;

            if t > 100.0 && beta < betaold && qfirst {
                ipz = 0;
                qfirst = false;
            } else if qfirst {
                ipz = 1;
            }

            betaold = beta;

            array_distrtkm.push(distrtkm);
            array_altkm.push(altkm);
            array_vkm.push(vkm);
        }
    }

    Results {
        distrtkm: array_distrtkm,
        altkm: array_altkm,
        vkm: array_vkm,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c15l8_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.distrtkm.clone(),
        results.altkm.clone(),
        results.vkm.clone(),
    ])?;

    let plot_file = format!("{}/c15l8_trajectory.png", output_dir);
    let config = PlotConfig::new("High Trajectory (Lob)")
        .with_labels("Downrange (km)", "Altitude (km)");
    let series = vec![
        Series::new(results.distrtkm.clone(), results.altkm.clone())
            .with_color(plotters::prelude::BLUE),
    ];
    line_plot(&plot_file, &config, &series).ok();

    println!("C15L8: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c15l8_runs() {
        let results = run();
        assert!(!results.distrtkm.is_empty());
    }
}
