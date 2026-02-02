//! Chapter 26, Lesson 4: Optimal Control with Gyro Dynamics
//!
//! Simulates optimal control with actuator and gyro dynamics.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub phi: Vec<f64>,     // Angle
    pub phid: Vec<f64>,    // Angular rate
    pub del: Vec<f64>,     // Deflection
}

/// Run the C26L4 simulation
pub fn run() -> Results {
    // QGYR=0; then QGYR=1
    let qgyr: i32 = 1;
    let wg: f64 = 200.0;
    let zg: f64 = 0.5;
    let wact: f64 = 100.0;
    let zact: f64 = 0.65;
    let xkd: f64 = 9000.0;
    let wr: f64 = 2.0;
    let c1: f64 = 3.0;
    let c2: f64 = 0.127;
    let c3: f64 = 8.81;
    let c4: f64 = 0.0309;

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
        let (_delc, deldd, phidd, phimddd) = if qgyr == 1 {
            let delc = -c1 * phim - c2 * phimd - c3 * del - c4 * deld;
            let deldd = wact * wact * (delc - del - 2.0 * zact * deld / wact);
            let phidd = xkd * del - wr * phid;
            let phimddd = wg * wg * (phid - phimd - 2.0 * zg * phimdd / wg);
            (delc, deldd, phidd, phimddd)
        } else {
            let delc = -c1 * phi - c2 * phid - c3 * del - c4 * deld;
            let deldd = wact * wact * (delc - del - 2.0 * zact * deld / wact);
            let phidd = xkd * del - wr * phid;
            let phimddd = 0.0;
            (delc, deldd, phidd, phimddd)
        };

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
        let (_delc, _deldd, phidd, _phimddd) = if qgyr == 1 {
            let delc = -c1 * phim - c2 * phimd - c3 * del - c4 * deld;
            let deldd = wact * wact * (delc - del - 2.0 * zact * deld / wact);
            let phidd = xkd * del - wr * phid;
            let phimddd = wg * wg * (phid - phimd - 2.0 * zg * phimdd / wg);
            (delc, deldd, phidd, phimddd)
        } else {
            let delc = -c1 * phi - c2 * phid - c3 * del - c4 * deld;
            let deldd = wact * wact * (delc - del - 2.0 * zact * deld / wact);
            let phidd = xkd * del - wr * phid;
            let phimddd = 0.0;
            (delc, deldd, phidd, phimddd)
        };

        // RK2 averaging
        phi = 0.5 * (phiold + phi + h * phid);
        phid = 0.5 * (phidold + phid + h * phidd);
        // Note: The MATLAB code does not RK2 average del, deld, phim, phimd, phimdd explicitly
        // but we can see from the structure it only averages phi and phid

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

    let data_file = format!("{}/c26l4_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.phi.clone(),
        results.phid.clone(),
        results.del.clone(),
    ])?;

    let plot_file = format!("{}/c26l4_phi.png", output_dir);
    let config = PlotConfig::new("Optimal Control with Gyro - PHI")
        .with_labels("Time (s)", "PHI (deg)");

    let series = vec![
        Series::new(results.time.clone(), results.phi.clone())
            .with_color(plotters::prelude::BLUE)
            .with_label("PHI"),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C26L4: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c26l4_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
