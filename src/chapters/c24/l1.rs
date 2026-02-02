//! Chapter 24, Lesson 1: Missile Aerodynamics
//!
//! Same aerodynamic analysis as C21L1 but for chapter 24 context.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub xnlg: Vec<f64>,
    pub alfdeg: Vec<f64>,
}

/// Run the C24L1 simulation
pub fn run() -> Results {
    let vm: f64 = 3000.0;
    let del: f64 = 5.0 / 57.3;
    let alt: f64 = 0.0;
    let a: f64 = 1000.0;
    let diam: f64 = 1.0;
    let fr: f64 = 3.0;
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

    let mut thd: f64 = 0.0;
    let mut alf: f64 = 0.0;
    let mut t: f64 = 0.0;
    let h: f64 = 0.0025;
    let mut s: f64 = 0.0;

    let mut array_t = Vec::new();
    let mut array_xnlg = Vec::new();
    let mut array_alfdeg = Vec::new();

    while t < 1.99999 {
        let thdold = thd;
        let alfold = alf;

        // First derivative evaluation
        let cn = 2.0 * alf + 1.5 * splan * alf * alf / sref
            + 8.0 * swing * alf / (b * sref)
            + 8.0 * stail * (alf + del) / (b * sref);
        let cm = 2.0 * alf * tmp4 + 1.5 * splan * alf * alf * tmp3 / sref
            + 8.0 * swing * alf * tmp1 / (b * sref)
            + 8.0 * stail * (alf + del) * tmp2 / (b * sref);

        let thdd = q * sref * diam * cm / xiyy;
        let xnl = 32.2 * q * sref * cn / wgt;
        let alfd = thd - xnl / vm;

        // Euler step
        thd += h * thdd;
        alf += h * alfd;
        t += h;

        // Second derivative for RK2
        let cn = 2.0 * alf + 1.5 * splan * alf * alf / sref
            + 8.0 * swing * alf / (b * sref)
            + 8.0 * stail * (alf + del) / (b * sref);
        let cm = 2.0 * alf * tmp4 + 1.5 * splan * alf * alf * tmp3 / sref
            + 8.0 * swing * alf * tmp1 / (b * sref)
            + 8.0 * stail * (alf + del) * tmp2 / (b * sref);

        let thdd = q * sref * diam * cm / xiyy;
        let xnl = 32.2 * q * sref * cn / wgt;
        let alfd = thd - xnl / vm;

        // RK2 averaging
        thd = 0.5 * (thdold + thd + h * thdd);
        alf = 0.5 * (alfold + alf + h * alfd);

        s += h;
        if s >= 0.0099999 {
            s = 0.0;
            array_t.push(t);
            array_xnlg.push(xnl / 32.2);
            array_alfdeg.push(alf * 57.3);
        }
    }

    Results {
        time: array_t,
        xnlg: array_xnlg,
        alfdeg: array_alfdeg,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c24l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.xnlg.clone(),
        results.alfdeg.clone(),
    ])?;

    let plot_file = format!("{}/c24l1_accel.png", output_dir);
    let config = PlotConfig::new("Missile Aerodynamics - Acceleration")
        .with_labels("Time (Sec)", "Missile Acceleration (G)");

    let series = vec![
        Series::new(results.time.clone(), results.xnlg.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C24L1: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c24l1_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
