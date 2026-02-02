//! Chapter 22, Lesson 1: Missile Trim Analysis
//!
//! Computes trim conditions and linearized transfer function response
//! for missile acceleration control.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub xnl: Vec<f64>,     // Actual acceleration (G)
    pub xncg: Vec<f64>,    // Commanded acceleration (G)
}

/// Run the C22L1 simulation
pub fn run() -> Results {
    let vm: f64 = 3000.0;
    let xncg: f64 = 10.0;  // Commanded acceleration (G)
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

    // Stability derivatives
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
    let xkdc = 1.0 / xk1;

    let mut e: f64 = 0.0;
    let mut ed: f64 = 0.0;
    let mut t: f64 = 0.0;
    let h: f64 = 0.0001;
    let mut s: f64 = 0.0;

    let mut array_t = Vec::new();
    let mut array_xnl = Vec::new();
    let mut array_xncg = Vec::new();

    while t < 1.99999 {
        let eold = e;
        let edold = ed;

        // First derivative evaluation
        let del = xkdc * xncg;
        let edd = waf * waf * (del - e - 2.0 * zaf * ed / waf);
        let _xnl = xk1 * (e - edd / (wz * wz));

        // Euler step
        e += h * ed;
        ed += h * edd;
        t += h;

        // Second derivative for RK2
        let del = xkdc * xncg;
        let edd = waf * waf * (del - e - 2.0 * zaf * ed / waf);
        let xnl = xk1 * (e - edd / (wz * wz));

        // RK2 averaging
        e = 0.5 * (eold + e + h * ed);
        ed = 0.5 * (edold + ed + h * edd);

        s += h;
        if s >= 0.0099999 {
            s = 0.0;
            array_t.push(t);
            array_xnl.push(xnl);
            array_xncg.push(xncg);
        }
    }

    Results {
        time: array_t,
        xnl: array_xnl,
        xncg: array_xncg,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c22l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.xnl.clone(),
        results.xncg.clone(),
    ])?;

    let plot_file = format!("{}/c22l1_accel.png", output_dir);
    let config = PlotConfig::new("Trim Analysis - Acceleration Response")
        .with_labels("Time (Sec)", "Missile Acceleration (G)");

    let series = vec![
        Series::new(results.time.clone(), results.xnl.clone())
            .with_color(plotters::prelude::BLUE)
            .with_label("Actual"),
        Series::new(results.time.clone(), results.xncg.clone())
            .with_color(plotters::prelude::RED)
            .with_label("Commanded"),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C22L1: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c22l1_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
