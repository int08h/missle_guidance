//! Chapter 26, Lesson 2: Optimal Control Simulation
//!
//! Simulates optimal control with actuator dynamics.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub phi: Vec<f64>,     // Angle
    pub phid: Vec<f64>,    // Angular rate
    pub del: Vec<f64>,     // Deflection
}

/// Run the C26L2 simulation
pub fn run() -> Results {
    // QACT=1; then QACT=0
    let qact: i32 = 0;
    let wact: f64 = 100.0;
    let zact: f64 = 0.65;
    let xkd: f64 = 9000.0;
    let wrr: f64 = 2.0;
    let c1: f64 = 3.0;
    let c2: f64 = 0.103;

    let mut t: f64 = 0.0;
    let h: f64 = 0.0002;
    let mut s: f64 = 0.0;
    let mut phi: f64 = 10.0;
    let mut phid: f64 = 0.0;
    let mut del: f64 = 0.0;
    let mut deld: f64 = 0.0;

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

        // First derivative evaluation
        let delc = -c1 * phi - c2 * phid;
        let deldd = wact * wact * (delc - del - 2.0 * zact * deld / wact);
        let e = if qact == 1 { del } else { delc };
        let phidd = xkd * e - wrr * phid;

        // Euler step
        phi += h * phid;
        phid += h * phidd;
        del += h * deld;
        deld += h * deldd;
        t += h;

        // Second derivative for RK2
        let delc = -c1 * phi - c2 * phid;
        let _deldd = wact * wact * (delc - del - 2.0 * zact * deld / wact);
        let e = if qact == 1 { del } else { delc };
        let phidd = xkd * e - wrr * phid;

        // RK2 averaging
        phi = 0.5 * (phiold + phi + h * phid);
        phid = 0.5 * (phidold + phid + h * phidd);
        // Note: del and deld are not RK2 averaged in MATLAB for this listing
        // They are just stepped forward but we record values at sample times

        // Sample output (S>=.0-9999 is essentially always true, so sample every step)
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

    let data_file = format!("{}/c26l2_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.phi.clone(),
        results.phid.clone(),
        results.del.clone(),
    ])?;

    let plot_file = format!("{}/c26l2_phi.png", output_dir);
    let config = PlotConfig::new("Optimal Control - PHI")
        .with_labels("Time (s)", "PHI (deg)");

    let series = vec![
        Series::new(results.time.clone(), results.phi.clone())
            .with_color(plotters::prelude::BLUE)
            .with_label("PHI"),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C26L2: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c26l2_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
