//! Chapter 15, Lesson 4: Ballistic Trajectory with Given Distance
//!
//! Ballistic trajectory to hit target at 6000 nmi with optimal launch angle.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub xnm: Vec<f64>,
    pub ynm: Vec<f64>,
    pub distnm: Vec<f64>,
    pub altnm: Vec<f64>,
}

/// Run the C15L4 simulation
pub fn run() -> Results {
    let h: f64 = 0.01;
    let a: f64 = 2.0926e7;
    let gm: f64 = 1.4077e16;
    let gamdeg: f64 = 23.0;
    let gam = gamdeg / 57.3;
    let distnm_target: f64 = 6000.0;
    let angdeg: f64 = 0.0;
    let ang = angdeg / 57.3;
    let phi = distnm_target * 6076.0 / a;
    let altnm_init: f64 = 0.0;
    let mut alt = altnm_init * 6076.0;
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
    let mut s: f64 = 0.0;

    let mut array_t = Vec::new();
    let mut array_xnm = Vec::new();
    let mut array_ynm = Vec::new();
    let mut array_distnm = Vec::new();
    let mut array_altnm = Vec::new();

    while alt >= 0.0 {
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

        s += h;
        if s >= 9.99999 {
            s = 0.0;
            let xnm = x / 6076.0;
            let ynm = y / 6076.0;
            alt = (x * x + y * y).sqrt() - a;
            let altnm = alt / 6076.0;
            let r = (x * x + y * y).sqrt();
            let rf = (xfirst * xfirst + yfirst * yfirst).sqrt();
            let cbeta = (x * xfirst + y * yfirst) / (r * rf);
            let beta = cbeta.acos();
            let distnm = a * beta / 6076.0;

            array_t.push(t);
            array_xnm.push(xnm);
            array_ynm.push(ynm);
            array_distnm.push(distnm);
            array_altnm.push(altnm);
        }
    }

    Results {
        time: array_t,
        xnm: array_xnm,
        ynm: array_ynm,
        distnm: array_distnm,
        altnm: array_altnm,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c15l4_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.xnm.clone(),
        results.ynm.clone(),
        results.distnm.clone(),
        results.altnm.clone(),
    ])?;

    let plot_file = format!("{}/c15l4_trajectory.png", output_dir);
    let config = PlotConfig::new("Ballistic Trajectory")
        .with_labels("Downrange (Nmi)", "Altitude (Nmi)");
    let series = vec![
        Series::new(results.distnm.clone(), results.altnm.clone())
            .with_color(plotters::prelude::BLUE),
    ];
    line_plot(&plot_file, &config, &series).ok();

    println!("C15L4: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c15l4_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
