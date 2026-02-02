//! Chapter 15, Lesson 3: Circular Orbit Simulation
//!
//! Satellite in circular orbit at 1000 nmi altitude.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub xnm: Vec<f64>,
    pub ynm: Vec<f64>,
}

/// Run the C15L3 simulation
pub fn run() -> Results {
    let h: f64 = 0.01;
    let a: f64 = 2.0926e7;
    let gm: f64 = 1.4077e16;
    let gam: f64 = 0.0;
    let altnm: f64 = 1000.0;
    let alt = altnm * 6076.0;
    let xlam: f64 = 1.0;
    let v = (gm * xlam / (a + alt)).sqrt();
    let angdeg: f64 = 90.0;
    let ang = angdeg / 57.3;

    let vrx = v * (1.5708 - gam / 57.3 + ang).cos();
    let vry = v * (1.5708 - gam / 57.3 + ang).sin();

    let mut x = (a + alt) * ang.cos();
    let mut y = (a + alt) * ang.sin();
    let mut x1 = vrx;
    let mut y1 = vry;
    let mut t: f64 = 0.0;
    let mut s: f64 = 0.0;
    let tf: f64 = 30000.0;

    let mut array_t = Vec::new();
    let mut array_xnm = Vec::new();
    let mut array_ynm = Vec::new();

    while t < tf {
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
        if s >= 49.99999 {
            s = 0.0;
            let xnm = x / 6076.0;
            let ynm = y / 6076.0;

            array_t.push(t);
            array_xnm.push(xnm);
            array_ynm.push(ynm);
        }
    }

    Results {
        time: array_t,
        xnm: array_xnm,
        ynm: array_ynm,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c15l3_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.xnm.clone(),
        results.ynm.clone(),
    ])?;

    let plot_file = format!("{}/c15l3_orbit.png", output_dir);
    let config = PlotConfig::new("Circular Orbit")
        .with_labels("X (Nmi)", "Y (Nmi)");
    let series = vec![
        Series::new(results.xnm.clone(), results.ynm.clone())
            .with_color(plotters::prelude::BLUE),
    ];
    line_plot(&plot_file, &config, &series).ok();

    println!("C15L3: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c15l3_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
