//! Chapter 15, Lesson 5: Minimum Energy Ballistic Trajectory
//!
//! Ballistic trajectory with minimum energy flight path angle.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub distrtkm: Vec<f64>,
    pub altkm: Vec<f64>,
    pub velkm: Vec<f64>,
}

/// Run the C15L5 simulation
pub fn run() -> Results {
    let qmin: i32 = 1;
    let h: f64 = 0.01;
    let a: f64 = 2.0926e7;
    let gm: f64 = 1.4077e16;
    let distkm: f64 = 10000.0;
    let angdeg: f64 = 0.0;
    let ang = angdeg / 57.3;
    let phi = distkm * 3280.0 / a;

    let gam = if qmin == 1 {
        std::f64::consts::PI / 4.0 - phi / 4.0
    } else {
        20.0 / 57.3
    };

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

    let mut array_t = Vec::new();
    let mut array_distrtkm = Vec::new();
    let mut array_altkm = Vec::new();
    let mut array_velkm = Vec::new();

    while !(t > 10.0 && alt < 0.0) {
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
        alt = (x * x + y * y).sqrt() - a;

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
        if scount >= 9.9999 {
            scount = 0.0;
            let altkm = alt / 3280.0;
            let r = (x * x + y * y).sqrt();
            let rf = (xfirst * xfirst + yfirst * yfirst).sqrt();
            let mut cbeta = (x * xfirst + y * yfirst) / (r * rf);
            if cbeta > 1.0 {
                cbeta = 1.0;
            }
            let beta = cbeta.acos();
            let distrtkm = a * beta / 3280.0;
            let velkm = (x1 * x1 + y1 * y1).sqrt() / 3280.0;

            array_t.push(t);
            array_distrtkm.push(distrtkm);
            array_altkm.push(altkm);
            array_velkm.push(velkm);
        }
    }

    Results {
        time: array_t,
        distrtkm: array_distrtkm,
        altkm: array_altkm,
        velkm: array_velkm,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c15l5_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.distrtkm.clone(),
        results.altkm.clone(),
        results.velkm.clone(),
    ])?;

    let plot_file = format!("{}/c15l5_trajectory.png", output_dir);
    let config = PlotConfig::new("Minimum Energy Trajectory")
        .with_labels("Downrange (km)", "Altitude (km)");
    let series = vec![
        Series::new(results.distrtkm.clone(), results.altkm.clone())
            .with_color(plotters::prelude::BLUE),
    ];
    line_plot(&plot_file, &config, &series).ok();

    let plot_file2 = format!("{}/c15l5_velocity.png", output_dir);
    let config2 = PlotConfig::new("Velocity vs Time")
        .with_labels("Time (s)", "Velocity (km/s)");
    let series2 = vec![
        Series::new(results.time.clone(), results.velkm.clone())
            .with_color(plotters::prelude::BLUE),
    ];
    line_plot(&plot_file2, &config2, &series2).ok();

    println!("C15L5: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c15l5_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
