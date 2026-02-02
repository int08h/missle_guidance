//! Chapter 29, Lesson 3: Normalized Miss vs X Parameter
//!
//! Computes maximum normalized miss distance as a function of the X parameter
//! (weave frequency or time constant ratio).

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub x: Vec<f64>,
    pub xmweavemax: Vec<f64>,
}

/// Run the C29L3 simulation
pub fn run() -> Results {
    let vc: f64 = 4000.0;
    let xnt: f64 = 32.2;
    let xnp: f64 = 3.0;
    let xnclim: f64 = 99999.0;
    let h: f64 = 0.01;

    let mut array_x = Vec::new();
    let mut array_xmweavemax = Vec::new();

    let mut x_param = 0.1;
    while x_param <= 4.0 + 1e-9 {
        let (w, tau) = if x_param < 0.5 {
            (1.0, x_param / 1.0)
        } else {
            (x_param, 1.0)
        };

        let mut xmweaveold: f64 = 0.0;
        let mut xmweavemax: f64 = 0.0;

        let mut tf = 0.2;
        while tf <= 20.0 + 1e-9 {
            let mut y: f64 = 0.0;
            let mut yd: f64 = 0.0;
            let mut xnl: f64 = 0.0;
            let mut d: f64 = 0.0;
            let mut elamdh: f64 = 0.0;
            let mut x4: f64 = 0.0;
            let mut x5: f64 = 0.0;
            let mut t: f64 = 0.0;

            while t <= tf - 1e-5 {
                let yold = y;
                let ydold = yd;
                let xnlold = xnl;
                let dold = d;
                let elamdhold = elamdh;
                let x4old = x4;
                let x5old = x5;

                // First derivative evaluation
                let ytdd = xnt * (w * t).sin();
                let tgo = tf - t + 0.00001;
                let xlam = y / (vc * tgo);
                let dd = 5.0 * (xlam - d) / tau;
                let elamdhd = 5.0 * (dd - elamdh) / tau;
                let mut xnc = xnp * vc * elamdh;
                if xnc > xnclim {
                    xnc = xnclim;
                }
                if xnc < -xnclim {
                    xnc = -xnclim;
                }
                let x4d = 5.0 * (xnc - x4) / tau;
                let x5d = 5.0 * (x4 - x5) / tau;
                let xnld = 5.0 * (x5 - xnl) / tau;
                let ydd = ytdd - xnl;

                // Euler step
                y += h * yd;
                yd += h * ydd;
                xnl += h * xnld;
                elamdh += h * elamdhd;
                d += h * dd;
                x4 += h * x4d;
                x5 += h * x5d;
                t += h;

                // Second derivative evaluation
                let ytdd = xnt * (w * t).sin();
                let tgo = tf - t + 0.00001;
                let xlam = y / (vc * tgo);
                let dd = 5.0 * (xlam - d) / tau;
                let elamdhd = 5.0 * (dd - elamdh) / tau;
                let mut xnc = xnp * vc * elamdh;
                if xnc > xnclim {
                    xnc = xnclim;
                }
                if xnc < -xnclim {
                    xnc = -xnclim;
                }
                let x4d = 5.0 * (xnc - x4) / tau;
                let x5d = 5.0 * (x4 - x5) / tau;
                let xnld = 5.0 * (x5 - xnl) / tau;
                let ydd = ytdd - xnl;

                // RK2 averaging
                y = 0.5 * (yold + y + h * yd);
                yd = 0.5 * (ydold + yd + h * ydd);
                xnl = 0.5 * (xnlold + xnl + h * xnld);
                d = 0.5 * (dold + d + h * dd);
                elamdh = 0.5 * (elamdhold + elamdh + h * elamdhd);
                x4 = 0.5 * (x4old + x4 + h * x4d);
                x5 = 0.5 * (x5old + x5 + h * x5d);
            }

            let xmweave = y;
            if xmweave > xmweaveold && xmweave > xmweavemax && tf > 10.0 {
                xmweavemax = xmweave;
            }
            xmweaveold = xmweave;

            tf += 0.2;
        }

        if x_param < 0.5 {
            xmweavemax /= tau * tau;
        }

        array_x.push(x_param);
        array_xmweavemax.push(xmweavemax);

        x_param += 0.1;
    }

    Results {
        x: array_x,
        xmweavemax: array_xmweavemax,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c29l3_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.x.clone(),
        results.xmweavemax.clone(),
    ])?;

    let plot_file = format!("{}/c29l3_miss.png", output_dir);
    let config = PlotConfig::new("Normalized Miss vs X Parameter")
        .with_labels("X", "Normalized Miss");

    let series = vec![
        Series::new(results.x.clone(), results.xmweavemax.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C29L3: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c29l3_runs() {
        let results = run();
        assert!(!results.x.is_empty());
    }
}
