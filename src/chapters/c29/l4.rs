//! Chapter 29, Lesson 4: Miss Distance vs Flight Time with Different Guidance Laws
//!
//! Computes miss distance as a function of flight time using different
//! guidance laws (PN, APN with target model).

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub tf: Vec<f64>,
    pub y: Vec<f64>,
}

/// Run the C29L4 simulation
pub fn run() -> Results {
    let vc: f64 = 4000.0;
    let xnt: f64 = 193.2;
    let xnp: f64 = 3.0;
    let xnclim: f64 = 99999999.0;
    let tau: f64 = 0.25;
    let w: f64 = 2.0;
    let wh: f64 = 2.0;
    let apn: i32 = 1;
    let h: f64 = 0.01;

    let mut array_tf = Vec::new();
    let mut array_y = Vec::new();

    let mut tf = 0.1;
    while tf <= 10.0 + 1e-9 {
        let mut y: f64 = 0.0;
        let mut yd: f64 = 0.0;
        let mut xnl: f64 = 0.0;
        let mut d: f64 = 0.0;
        let mut elamdh: f64 = 0.0;
        let mut x4: f64 = 0.0;
        let mut x5: f64 = 0.0;
        let mut t: f64 = 0.0;

        while t <= tf - 0.0001 {
            let yold = y;
            let ydold = yd;
            let xnlold = xnl;
            let dold = d;
            let elamdhold = elamdh;
            let x4old = x4;
            let x5old = x5;

            // First derivative evaluation
            let ytdd = xnt * (w * t).sin();
            let ytddd = w * xnt * (w * t).cos();
            let tgo = tf - t + 0.00001;
            let xlam = y / (vc * tgo);
            let dd = 5.0 * (xlam - d) / tau;
            let elamdhd = 5.0 * (dd - elamdh) / tau;

            let xnc = if apn == 1 {
                xnp * vc * elamdh
            } else if apn == 2 {
                let xp = wh * tgo;
                xnp * vc * elamdh + xnp * ytdd * (1.0 - xp.cos()) / (xp * xp)
                    + xnp * ytddd * (xp - xp.sin()) / (xp * xp * wh)
            } else {
                let xs = tgo / tau;
                let xp = wh * tgo;
                let top = 6.0 * xs * xs * ((-xs).exp() - 1.0 + xs);
                let bot1 = 2.0 * xs.powi(3) + 3.0 + 6.0 * xs - 6.0 * xs * xs;
                let bot2 = -12.0 * xs * (-xs).exp() - 3.0 * (-2.0 * xs).exp();
                let xnpp = top / (0.0001 + bot1 + bot2);
                xnpp * vc * elamdh + xnpp * ytdd * (1.0 - xp.cos()) / (xp * xp)
                    + xnpp * ytddd * (xp - xp.sin()) / (xp * xp * wh)
                    - xnpp * xnl * tau * tau * ((-xs).exp() + xs - 1.0) / (tgo * tgo)
            };

            let xnc = xnc.clamp(-xnclim, xnclim);

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
            let ytddd = w * xnt * (w * t).cos();
            let tgo = tf - t + 0.00001;
            let xlam = y / (vc * tgo);
            let dd = 5.0 * (xlam - d) / tau;
            let elamdhd = 5.0 * (dd - elamdh) / tau;

            let xnc = if apn == 1 {
                xnp * vc * elamdh
            } else if apn == 2 {
                let xp = wh * tgo;
                xnp * vc * elamdh + xnp * ytdd * (1.0 - xp.cos()) / (xp * xp)
                    + xnp * ytddd * (xp - xp.sin()) / (xp * xp * wh)
            } else {
                let xs = tgo / tau;
                let xp = wh * tgo;
                let top = 6.0 * xs * xs * ((-xs).exp() - 1.0 + xs);
                let bot1 = 2.0 * xs.powi(3) + 3.0 + 6.0 * xs - 6.0 * xs * xs;
                let bot2 = -12.0 * xs * (-xs).exp() - 3.0 * (-2.0 * xs).exp();
                let xnpp = top / (0.0001 + bot1 + bot2);
                xnpp * vc * elamdh + xnpp * ytdd * (1.0 - xp.cos()) / (xp * xp)
                    + xnpp * ytddd * (xp - xp.sin()) / (xp * xp * wh)
                    - xnpp * xnl * tau * tau * ((-xs).exp() + xs - 1.0) / (tgo * tgo)
            };

            let xnc = xnc.clamp(-xnclim, xnclim);

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

        array_tf.push(tf);
        array_y.push(y);

        tf += 0.1;
    }

    Results {
        tf: array_tf,
        y: array_y,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c29l4_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.tf.clone(),
        results.y.clone(),
    ])?;

    let plot_file = format!("{}/c29l4_miss.png", output_dir);
    let config = PlotConfig::new("Miss Distance vs Flight Time")
        .with_labels("Flight Time (Sec)", "Miss (Ft)");

    let series = vec![
        Series::new(results.tf.clone(), results.y.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C29L4: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c29l4_runs() {
        let results = run();
        assert!(!results.tf.is_empty());
    }
}
