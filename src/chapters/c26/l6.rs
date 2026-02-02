//! Chapter 26, Lesson 6: Optimal Control with Gyro Dynamics - Extended
//!
//! Simulates optimal control with actuator and gyro dynamics using 6 gain terms.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub phi: Vec<f64>,     // Angle
    pub phid: Vec<f64>,    // Angular rate
    pub del: Vec<f64>,     // Deflection
}

/// Run the C26L6 simulation
pub fn run() -> Results {
    let wg: f64 = 200.0;
    let zg: f64 = 0.5;
    let wact: f64 = 100.0;
    let zact: f64 = 0.65;
    let xkd: f64 = 9000.0;
    let wr: f64 = 2.0;
    let c1: f64 = 3.0;
    let c2: f64 = 0.127;
    let c3: f64 = 8.84;
    let c4: f64 = 0.031;
    let c5: f64 = 0.0152;
    let c6: f64 = 0.0000768;

    let mut t: f64 = 0.0;
    let h: f64 = 0.0002;
    let mut s: f64 = 0.0;
    let mut phi: f64 = 10.0;
    let mut phid: f64 = 0.0;
    let mut del: f64 = 0.0;
    let mut deld: f64 = 0.0;
    let mut phim: f64 = 10.0;
    let mut phimd: f64 = 0.0;
    let mut phimdd: f64 = 0.0;

    let mut array_t = Vec::new();
    let mut array_phi = Vec::new();
    let mut array_phid = Vec::new();
    let mut array_del = Vec::new();

    while t <= 0.5 {
        s += h;
        let phiold = phi;
        let phidold = phid;
        let _delold = del;
        let _deldold = deld;
        let _phimold = phim;
        let _phimdold = phimd;
        let _phimddold = phimdd;

        // First derivative evaluation
        let delc = -c1 * phi - c2 * phid - c3 * del - c4 * deld - c5 * phimd - c6 * phimdd;
        let deldd = wact * wact * (delc - del - 2.0 * zact * deld / wact);
        let e = del;
        let phidd = xkd * e - wr * phid;
        let phimddd = wg * wg * (phid - phimd - 2.0 * zg * phimdd / wg);

        // Euler step
        phi += h * phid;
        phid += h * phidd;
        del += h * deld;
        deld += h * deldd;
        phim += h * phimd;
        phimd += h * phimdd;
        phimdd += h * phimddd;
        t += h;

        // Second derivative for RK2
        let delc = -c1 * phi - c2 * phid - c3 * del - c4 * deld - c5 * phimd - c6 * phimdd;
        let _deldd = wact * wact * (delc - del - 2.0 * zact * deld / wact);
        let e = del;
        let phidd = xkd * e - wr * phid;
        let _phimddd = wg * wg * (phid - phimd - 2.0 * zg * phimdd / wg);

        // RK2 averaging
        phi = 0.5 * (phiold + phi + h * phid);
        phid = 0.5 * (phidold + phid + h * phidd);

        // Sample output
        if s >= 0.0 - 9999.0 {
            s = 0.0;
            array_t.push(t);
            array_phi.push(phi);
            array_phid.push(phid);
            array_del.push(del);
        }
    }

    Results {
        time: array_t,
        phi: array_phi,
        phid: array_phid,
        del: array_del,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c26l6_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.phi.clone(),
        results.phid.clone(),
        results.del.clone(),
    ])?;

    let plot_file = format!("{}/c26l6_phi.png", output_dir);
    let config = PlotConfig::new("Optimal Control with Gyro (6 gains) - PHI")
        .with_labels("Time (s)", "PHI (deg)");

    let series = vec![
        Series::new(results.time.clone(), results.phi.clone())
            .with_color(plotters::prelude::BLUE)
            .with_label("PHI"),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C26L6: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c26l6_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
