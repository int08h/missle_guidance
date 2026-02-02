//! Chapter 22, Lesson 5: Numerical Frequency Response
//!
//! Computes frequency response (gain and phase) using time-domain
//! simulation rather than analytical transfer function.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub freq: Vec<f64>,    // Frequency (rad/sec)
    pub gain: Vec<f64>,    // Gain (dB)
    pub phase: Vec<f64>,   // Phase (deg)
}

/// Run the C22L5 simulation
pub fn run() -> Results {
    let zact: f64 = 0.7;
    let wact: f64 = 150.0;
    let k3: f64 = -1.89;
    let ta: f64 = 0.457;
    let zaf: f64 = 0.058;
    let waf: f64 = 25.3;
    let kr: f64 = 0.1;
    let pi: f64 = 3.1416;
    let h: f64 = 0.0001;

    let mut array_w = Vec::new();
    let mut array_gain = Vec::new();
    let mut array_phase = Vec::new();

    for i in 2..=160 {
        let w = 10.0_f64.powf(0.025 * (i as f64) - 1.0);
        let period = 2.0 * pi / w;

        let mut t: f64 = 0.0;
        let mut s: f64 = 0.0;
        let mut e: f64 = 0.0;
        let mut ed: f64 = 0.0;
        let mut del: f64 = 0.0;
        let mut deld: f64 = 0.0;
        let mut p: f64 = 0.0;
        let mut q: f64 = 0.0;
        let mut pprev: f64 = 0.0;
        let mut qprev: f64 = 0.0;
        let mut delp: f64 = 0.0;
        let mut delq: f64 = 0.0;
        let mut delpold: f64 = 0.0;
        let mut delqold: f64 = 0.0;
        let mut deldelp: f64 = 100.0;

        // Run until convergence
        while !(t > 20.0 && deldelp.abs() < 0.0001) {
            let eold = e;
            let edold = ed;
            let delold = del;
            let deldold = deld;
            let pold = p;
            let qold = q;

            // First derivative evaluation
            let x = -(w * t).sin();
            let deldd = wact * wact * (x - del - 2.0 * zact * deld / wact);
            let edd = waf * waf * (del - e - 2.0 * zaf * ed / waf);
            let y = kr * k3 * (e + ta * ed);
            let pd = y * (w * t).sin();
            let qd = y * (w * t).cos();

            // Euler step
            e += h * ed;
            ed += h * edd;
            del += h * deld;
            deld += h * deldd;
            p += h * pd;
            q += h * qd;
            t += h;

            // Second derivative for RK2
            let x = -(w * t).sin();
            let deldd = wact * wact * (x - del - 2.0 * zact * deld / wact);
            let edd = waf * waf * (del - e - 2.0 * zaf * ed / waf);
            let y = kr * k3 * (e + ta * ed);
            let pd = y * (w * t).sin();
            let qd = y * (w * t).cos();

            // RK2 averaging
            e = 0.5 * (eold + e + h * ed);
            ed = 0.5 * (edold + ed + h * edd);
            del = 0.5 * (delold + del + h * deld);
            deld = 0.5 * (deldold + deld + h * deldd);
            p = 0.5 * (pold + p + h * pd);
            q = 0.5 * (qold + q + h * qd);

            s += h;
            if s >= period - 0.0001 {
                s = 0.0;
                delp = p - pprev;
                delq = q - qprev;
                pprev = p;
                qprev = q;
                deldelp = delpold - delp;
                let _deldelq = delqold - delq;
                delpold = delp;
                delqold = delq;
            }
        }

        let mut phase = 57.3 * delq.atan2(delp);
        if phase > 90.0 {
            phase -= 360.0;
        }

        let gain = 10.0 * ((delp * delp + delq * delq) * w * w / (pi * pi)).log10();

        array_w.push(w);
        array_gain.push(gain);
        array_phase.push(phase);
    }

    Results {
        freq: array_w,
        gain: array_gain,
        phase: array_phase,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c22l5_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.freq.clone(),
        results.gain.clone(),
        results.phase.clone(),
    ])?;

    let plot_file = format!("{}/c22l5_gain.png", output_dir);
    let config = PlotConfig::new("Numerical Frequency Response - Gain")
        .with_labels("Frequency (r/s)", "Gain (db)");

    let series = vec![
        Series::new(results.freq.clone(), results.gain.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    let plot_file2 = format!("{}/c22l5_phase.png", output_dir);
    let config2 = PlotConfig::new("Numerical Frequency Response - Phase")
        .with_labels("Frequency (r/s)", "Phase (deg)");

    let series2 = vec![
        Series::new(results.freq.clone(), results.phase.clone())
            .with_color(plotters::prelude::RED),
    ];

    line_plot(&plot_file2, &config2, &series2).ok();

    println!("C22L5: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c22l5_runs() {
        let results = run();
        assert!(!results.freq.is_empty());
        assert_eq!(results.freq.len(), 159);
    }
}
