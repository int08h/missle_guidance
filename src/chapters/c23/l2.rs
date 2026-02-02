//! Chapter 23, Lesson 2: Three-Loop Autopilot Frequency Response
//!
//! Bode plot analysis (gain and phase vs frequency) for the three-loop autopilot.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub w: Vec<f64>,
    pub gain: Vec<f64>,
    pub phase: Vec<f64>,
}

/// Run the C23L2 simulation
pub fn run() -> Results {
    // Missile parameters
    let fr: f64 = 3.0;
    let diam: f64 = 1.0;
    let xl: f64 = 20.0;
    let ctw: f64 = 0.0;
    let crw: f64 = 6.0;
    let hw: f64 = 2.0;
    let ctt: f64 = 0.0;
    let crt: f64 = 2.0;
    let ht: f64 = 2.0;
    let xn: f64 = 4.0;
    let xcg: f64 = 10.0;
    let xhl: f64 = 19.5;
    let wact: f64 = 150.0;
    let zact: f64 = 0.7;
    let vm: f64 = 3000.0;
    let xncg: f64 = 10.0;
    let wcr: f64 = 50.0;
    let zeta: f64 = 0.7;
    let tau: f64 = 0.3;
    let alt: f64 = 0.0;
    let a: f64 = 1000.0;
    let wgt: f64 = 1000.0;

    let rho = if alt <= 30000.0 {
        0.002378 * (-alt / 30000.0).exp()
    } else {
        0.0034 * (-alt / 22000.0).exp()
    };

    let swing = 0.5 * hw * (ctw + crw);
    let stail = 0.5 * ht * (ctt + crt);
    let sref = 3.1416 * diam * diam / 4.0;
    let xlp = fr * diam;
    let splan = (xl - xlp) * diam + 1.33 * xlp * diam / 2.0;
    let xcpn = 2.0 * xlp / 3.0;
    let an = 0.67 * xlp * diam;
    let ab = (xl - xlp) * diam;
    let xcpb = (0.67 * an * xlp + ab * (xlp + 0.5 * (xl - xlp))) / (an + ab);
    let xcpw = xlp + xn + 0.7 * crw - 0.2 * ctw;
    let xmach = vm / a;
    let xiyy = wgt * (3.0 * ((diam / 2.0).powi(2)) + xl * xl) / (12.0 * 32.2);

    let tmp1 = (xcg - xcpw) / diam;
    let tmp2 = (xcg - xhl) / diam;
    let tmp3 = (xcg - xcpb) / diam;
    let tmp4 = (xcg - xcpn) / diam;

    let b = (xmach * xmach - 1.0).sqrt();
    let q = 0.5 * rho * vm * vm;

    // Trim calculations
    let p1 = wgt * xncg / (q * sref);
    let y1 = 2.0 + 8.0 * swing / (b * sref) + 8.0 * stail / (b * sref);
    let y2 = 1.5 * splan / sref;
    let y3 = 8.0 * stail / (b * sref);
    let y4 = 2.0 * tmp4 + 8.0 * swing * tmp1 / (b * sref) + 8.0 * stail * tmp2 / (b * sref);
    let y5 = 1.5 * splan * tmp3 / sref;
    let y6 = 8.0 * stail * tmp2 / (b * sref);

    let p2 = y2 - y3 * y5 / y6;
    let p3 = y1 - y3 * y4 / y6;
    let alftr = (-p3 + (p3 * p3 + 4.0 * p2 * p1).sqrt()) / (2.0 * p2);

    let cna = 2.0 + 1.5 * splan * alftr / sref + 8.0 * swing / (b * sref) + 8.0 * stail / (b * sref);
    let cnd = 8.0 * stail / (b * sref);
    let cmap = 2.0 * tmp4 + 1.5 * splan * alftr * tmp3 / sref + 8.0 * swing * tmp1 / (b * sref);
    let cma = cmap + 8.0 * stail * tmp2 / (b * sref);
    let cmd = 8.0 * stail * tmp2 / (b * sref);

    let xma = q * sref * diam * cma / xiyy;
    let xmd = q * sref * diam * cmd / xiyy;
    let za = -32.2 * q * sref * cna / (wgt * vm);
    let zd = -32.2 * q * sref * cnd / (wgt * vm);

    let wz = ((xma * zd - za * xmd) / zd).sqrt();
    let waf = (-xma).sqrt();
    let zaf = 0.5 * waf * za / xma;
    let xk1 = -vm * (xma * zd - xmd * za) / (1845.0 * xma);
    let _xk2 = xk1;
    let ta = xmd / (xma * zd - xmd * za);
    let xk3 = 1845.0 * xk1 / vm;

    // Autopilot gains
    let w_design = (tau * wcr * (1.0 + 2.0 * zaf * waf / wcr) - 1.0) / (2.0 * zeta * tau);
    let w0 = w_design / (tau * wcr).sqrt();
    let z0 = 0.5 * w0 * (2.0 * zeta / w_design + tau - waf * waf / (w0 * w0 * wcr));
    let xkc = (-w0 * w0 / (wz * wz) - 1.0 + 2.0 * z0 * w0 * ta) / (1.0 - 2.0 * z0 * w0 * ta + w0 * w0 * ta * ta);
    let _xka = xk3 / (xk1 * xkc);
    let xk0 = -w_design * w_design / (tau * waf * waf);
    let _xk = xk0 / (xk1 * (1.0 + xkc));
    let _wi = xkc * ta * w0 * w0 / (1.0 + xkc + w0 * w0 / (wz * wz));
    let _xkr = _xk / (_xka * _wi);
    let _xkdc = 1.0 + 1845.0 / (_xka * vm);

    // Frequency response analysis
    let mut array_w = Vec::new();
    let mut array_gain = Vec::new();
    let mut array_phase = Vec::new();

    for i in 2..=160 {
        let w = 10.0_f64.powf(0.025 * (i as f64) - 1.0);

        let xmagtop = -xk0 * ((1.0 - (w / w0).powi(2)).powi(2) + (2.0 * z0 * w / w0).powi(2)).sqrt();
        let xmagbot = w * ((1.0 - (w / waf).powi(2)).powi(2) + (2.0 * zaf * w / waf).powi(2)).sqrt();
        let xmag = xmagtop / xmagbot;
        let xmagact = 1.0 / ((1.0 - w * w / (wact * wact)).powi(2) + (2.0 * zact * w / wact).powi(2)).sqrt();

        let phasetop = (2.0 * z0 * w / w0).atan2(1.0 - (w / w0).powi(2));
        let phasebot = (2.0 * zaf * w / waf).atan2(1.0 - (w / waf).powi(2));
        let phaseact = (2.0 * zact * w / wact).atan2(1.0 - w * w / (wact * wact));

        let gain = 20.0 * (xmag * xmagact).log10();
        let phase = -90.0 + 57.3 * (phasetop - phasebot - phaseact);

        array_w.push(w);
        array_gain.push(gain);
        array_phase.push(phase);
    }

    Results {
        w: array_w,
        gain: array_gain,
        phase: array_phase,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c23l2_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.w.clone(),
        results.gain.clone(),
        results.phase.clone(),
    ])?;

    // Gain plot
    let plot_file = format!("{}/c23l2_gain.png", output_dir);
    let config = PlotConfig::new("Three-Loop Autopilot - Frequency Response (Gain)")
        .with_labels("Frequency (Rad/Sec)", "Gain (dB)");

    let series = vec![
        Series::new(results.w.clone(), results.gain.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    // Phase plot
    let plot_file2 = format!("{}/c23l2_phase.png", output_dir);
    let config2 = PlotConfig::new("Three-Loop Autopilot - Frequency Response (Phase)")
        .with_labels("Frequency (Rad/Sec)", "Phase (Deg)");

    let series2 = vec![
        Series::new(results.w.clone(), results.phase.clone())
            .with_color(plotters::prelude::RED),
    ];

    line_plot(&plot_file2, &config2, &series2).ok();

    println!("C23L2: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c23l2_runs() {
        let results = run();
        assert_eq!(results.w.len(), 159);
    }
}
