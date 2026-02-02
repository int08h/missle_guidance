//! Chapter 38, Lesson 1: Trapezoidal Weave
//!
//! Models target acceleration with a trapezoidal weave pattern
//! and its Fourier series approximation.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

// Use MATLAB's PI value (3.14159) instead of std::f64::consts::PI
// to ensure numerical compatibility with the original code.
const PI: f64 = 3.14159;

pub struct Results {
    pub time: Vec<f64>,
    pub ytdd: Vec<f64>,
    pub xnth: Vec<f64>,
}

/// Run the C38L1 simulation
pub fn run() -> Results {
    let pz: f64 = 3.0;
    let tr: f64 = 0.5;
    let xl = pz / 2.0;
    let w = 2.0 * PI / pz;
    let xnt: f64 = 322.0;
    let x = pz / 2.0 - tr;
    let alf = PI * tr / (2.0 * xl);

    let mut array_t = Vec::new();
    let mut array_ytdd = Vec::new();
    let mut array_xnth = Vec::new();

    let mut t: f64 = 0.0;
    while t <= 6.0 {
        let tmod = t % pz;
        let ytdd = if tmod < tr / 2.0 {
            2.0 * xnt * tmod / tr
        } else if tmod < tr / 2.0 + x {
            xnt
        } else if tmod < 3.0 * tr / 2.0 + x {
            -2.0 * xnt * tmod / tr + 2.0 * xnt + 2.0 * xnt * x / tr
        } else if tmod < 3.0 * tr / 2.0 + 2.0 * x {
            -xnt
        } else {
            2.0 * xnt * tmod / tr - 4.0 * xnt - 4.0 * xnt * x / tr
        };

        // Fourier series approximation
        let xnth = 4.0 * xnt * (alf.sin() * (w * t).sin() + (3.0 * alf).sin() * (3.0 * w * t).sin() / 9.0) / (PI * alf);

        array_t.push(t);
        array_ytdd.push(ytdd);
        array_xnth.push(xnth);

        t += 0.1;
    }

    Results {
        time: array_t,
        ytdd: array_ytdd,
        xnth: array_xnth,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c38l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.ytdd.clone(),
        results.xnth.clone(),
    ])?;

    let plot_file = format!("{}/c38l1_weave.png", output_dir);
    let config = PlotConfig::new("Trapezoidal Weave")
        .with_labels("Time (Sec)", "Acceleration (ft/s^2)");

    let series = vec![
        Series::new(results.time.clone(), results.ytdd.clone())
            .with_color(plotters::prelude::BLUE)
            .with_label("Actual"),
        Series::new(results.time.clone(), results.xnth.clone())
            .with_color(plotters::prelude::RED)
            .with_label("Estimate"),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C38L1: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c38l1_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
