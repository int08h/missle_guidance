//! Chapter 4, Lesson 5: Shaping Filter Monte Carlo Results
//!
//! Monte Carlo simulation to determine miss distance statistics
//! using proportional navigation with a shaping filter.

use crate::plotting::{PlotConfig, scatter_plot};
use crate::save_data;
use rand::Rng;

/// Simulation results
pub struct Results {
    pub tf: Vec<f64>,
    pub sigma: Vec<f64>,
    pub xmean: Vec<f64>,
}

/// Run the C4L5 simulation
pub fn run() -> Results {
    let vc = 4000.0;
    let xnp = 3.0;
    let tau = 1.0;
    let num_runs = 50;
    let h = 0.01;

    let mut rng = rand::thread_rng();

    let mut array_tf = Vec::new();
    let mut array_sigma = Vec::new();
    let mut array_xmean = Vec::new();

    for tf in 1..=10 {
        let tf_f64 = tf as f64;
        let mut z = vec![0.0; num_runs];
        let mut z1 = 0.0;

        for z_item in z.iter_mut() {
            let tstart = tf_f64 * rng.gen::<f64>();
            let pz: f64 = rng.gen::<f64>() - 0.5;
            let coef = if pz > 0.0 { 1.0 } else { -1.0 };

            let mut y = 0.0;
            let mut yd = 0.0;
            let mut t = 0.0;
            let mut xnl = 0.0;

            while t <= (tf_f64 - 1e-5) {
                let xnt = if t < tstart { 0.0 } else { coef * 96.6 };

                let y_old = y;
                let yd_old = yd;
                let xnl_old = xnl;

                // First step of RK2
                let tgo = tf_f64 - t + 0.00001;
                let rtm = vc * tgo;
                let xlamd = (rtm * yd + y * vc) / (rtm * rtm);
                let xnc = xnp * vc * xlamd;
                let xnld = (xnc - xnl) / tau;
                let ydd = xnt - xnl;

                y += h * yd;
                yd += h * ydd;
                xnl += h * xnld;
                t += h;

                // Second step of RK2
                let tgo = tf_f64 - t + 0.00001;
                let rtm = vc * tgo;
                let xlamd = (rtm * yd + y * vc) / (rtm * rtm);
                let xnc = xnp * vc * xlamd;
                let xnld = (xnc - xnl) / tau;
                let ydd = xnt - xnl;

                y = 0.5 * (y_old + y + h * yd);
                yd = 0.5 * (yd_old + yd + h * ydd);
                xnl = 0.5 * (xnl_old + xnl + h * xnld);
            }

            *z_item = y;
            z1 += *z_item;
        }

        let xmean = z1 / num_runs as f64;

        // Calculate standard deviation
        let mut z1 = 0.0;
        for z_val in z.iter() {
            z1 += (*z_val - xmean).powi(2);
        }
        let sigma = if num_runs > 1 {
            (z1 / (num_runs - 1) as f64).sqrt()
        } else {
            0.0
        };

        array_tf.push(tf_f64);
        array_sigma.push(sigma);
        array_xmean.push(xmean);
    }

    Results {
        tf: array_tf,
        sigma: array_sigma,
        xmean: array_xmean,
    }
}

/// Run and save results to file
pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    // Save data file
    let data_file = format!("{}/c4l5_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.tf.clone(),
        results.sigma.clone(),
        results.xmean.clone(),
    ])?;

    // Create plot
    let plot_file = format!("{}/c4l5_plot.png", output_dir);
    let config = PlotConfig::new("Shaping filter Monte Carlo results")
        .with_labels("Time", "Standard Deviation / Mean")
        .with_y_range(0.0, 30.0);

    scatter_plot(&plot_file, &config, &results.tf, &results.sigma).ok();

    println!("C4L5: Simulation finished");
    println!("  Data saved to: {}", data_file);
    println!("  Plot saved to: {}", plot_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c4l5_runs() {
        let results = run();
        assert_eq!(results.tf.len(), 10);
        assert_eq!(results.sigma.len(), 10);
        assert_eq!(results.xmean.len(), 10);
    }
}
