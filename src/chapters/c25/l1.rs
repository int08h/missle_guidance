//! Chapter 25, Lesson 1: Flexible Body Effects
//!
//! Simulates rigid body and flexible body rates with actuator dynamics.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub thdrb: Vec<f64>,   // Rigid body rate
    pub thdfb: Vec<f64>,   // Flexible body rate
    pub thdtot: Vec<f64>,  // Total body rate
}

/// Run the C25L1 simulation
pub fn run() -> Results {
    let ts: f64 = 0.001;
    let xk3: f64 = -0.362;
    let ta: f64 = 0.831;
    let b11: f64 = 0.00461;
    let b12: f64 = 0.00136;
    let xkfb: f64 = 0.00134;
    let wzfb: f64 = 395.0;
    let z1: f64 = 0.015;
    let w1: f64 = 259.0;
    let tf: f64 = 1.0;
    let delc: f64 = 1.0;
    let wact: f64 = 100.0;
    let zact: f64 = 0.7;

    let mut t: f64 = 0.0;
    let mut s: f64 = 0.0;
    let h: f64 = 0.00001;
    let mut e: f64 = 0.0;
    let mut ed: f64 = 0.0;
    let mut e1: f64 = 0.0;
    let mut e1d: f64 = 0.0;
    let mut e1dd: f64 = 0.0;
    let mut del: f64 = 0.0;
    let mut deld: f64 = 0.0;

    let mut array_t = Vec::new();
    let mut array_thdrb = Vec::new();
    let mut array_thdfb = Vec::new();
    let mut array_thdtot = Vec::new();

    while t <= tf {
        let eold = e;
        let edold = ed;
        let e1old = e1;
        let e1dold = e1d;
        let e1ddold = e1dd;
        let delold = del;
        let deldold = deld;

        // First derivative evaluation
        let deldd = wact * wact * (delc - del - 2.0 * zact * deld / wact);
        let edd = (del - e - b11 * ed) / b12;
        let e1ddd = w1 * w1 * (deld - e1d - 2.0 * z1 * e1dd / w1);
        let thdrb = xk3 * (e + ta * ed);
        let thdfb = xkfb * (e1d + e1ddd / (wzfb * wzfb));
        let _thdtot = thdrb + thdfb;

        // Euler step
        e += h * ed;
        ed += h * edd;
        e1 += h * e1d;
        e1d += h * e1dd;
        e1dd += h * e1ddd;
        del += h * deld;
        deld += h * deldd;
        t += h;

        // Second derivative for RK2
        let deldd = wact * wact * (delc - del - 2.0 * zact * deld / wact);
        let edd = (del - e - b11 * ed) / b12;
        let e1ddd = w1 * w1 * (deld - e1d - 2.0 * z1 * e1dd / w1);
        let thdrb = xk3 * (e + ta * ed);
        let thdfb = xkfb * (e1d + e1ddd / (wzfb * wzfb));
        let thdtot = thdrb + thdfb;

        // RK2 averaging
        e = 0.5 * (eold + e + h * ed);
        ed = 0.5 * (edold + ed + h * edd);
        e1 = 0.5 * (e1old + e1 + h * e1d);
        e1d = 0.5 * (e1dold + e1d + h * e1dd);
        e1dd = 0.5 * (e1ddold + e1dd + h * e1ddd);
        del = 0.5 * (delold + del + h * deld);
        deld = 0.5 * (deldold + deld + h * deldd);

        s += h;
        if s >= ts - 0.00001 {
            s = 0.0;
            array_t.push(t);
            array_thdrb.push(thdrb);
            array_thdfb.push(thdfb);
            array_thdtot.push(thdtot);
        }
    }

    Results {
        time: array_t,
        thdrb: array_thdrb,
        thdfb: array_thdfb,
        thdtot: array_thdtot,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c25l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.thdrb.clone(),
        results.thdfb.clone(),
        results.thdtot.clone(),
    ])?;

    let plot_file = format!("{}/c25l1_rates.png", output_dir);
    let config = PlotConfig::new("Flexible Body Effects - Body Rates")
        .with_labels("Time (Sec)", "Body Rate (deg/s)");

    let series = vec![
        Series::new(results.time.clone(), results.thdrb.clone())
            .with_color(plotters::prelude::BLUE)
            .with_label("Rigid Body"),
        Series::new(results.time.clone(), results.thdfb.clone())
            .with_color(plotters::prelude::RED)
            .with_label("Flexible Body"),
        Series::new(results.time.clone(), results.thdtot.clone())
            .with_color(plotters::prelude::GREEN)
            .with_label("Total"),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C25L1: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c25l1_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
