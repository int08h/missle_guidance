//! Chapter 23, Lesson 3: Radome Slope Analysis
//!
//! Analysis of miss distance as a function of radome slope.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub r: Vec<f64>,
    pub xmfn: Vec<f64>,
    pub xmrn: Vec<f64>,
    pub xmgl: Vec<f64>,
    pub xmudnt: Vec<f64>,
    pub rms: Vec<f64>,
}

/// Run the C23L3 simulation
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
    let tf: f64 = 5.0;
    let vm: f64 = 3000.0;
    let xncg: f64 = 10.0;
    let wcr: f64 = 50.0;
    let zeta: f64 = 0.7;
    let tau: f64 = 0.3;
    let alt: f64 = 0.0;
    let xnt: f64 = 64.4;
    let xnp: f64 = 3.0;
    let vc: f64 = 4000.0;
    let t1: f64 = 0.1;
    let t2: f64 = 0.15;
    let phirn: f64 = 0.000002;
    let phigl: f64 = 20.0;
    let phifn: f64 = 0.00000008;
    let ra: f64 = 30000.0;
    let miss: i32 = 1;
    let tint: f64 = 0.0;
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
    let w = (tau * wcr * (1.0 + 2.0 * zaf * waf / wcr) - 1.0) / (2.0 * zeta * tau);
    let w0 = w / (tau * wcr).sqrt();
    let z0 = 0.5 * w0 * (2.0 * zeta / w + tau - waf * waf / (w0 * w0 * wcr));
    let xkc = (-w0 * w0 / (wz * wz) - 1.0 + 2.0 * z0 * w0 * ta) / (1.0 - 2.0 * z0 * w0 * ta + w0 * w0 * ta * ta);
    let xka = xk3 / (xk1 * xkc);
    let xk0 = -w * w / (tau * waf * waf);
    let _xk = xk0 / (xk1 * (1.0 + xkc));
    let wi = xkc * ta * w0 * w0 / (1.0 + xkc + w0 * w0 / (wz * wz));
    let xkr = _xk / (xka * wi);
    let xkdc = 1.0 + 1845.0 / (xka * vm);

    let h: f64 = 0.0001;

    let mut array_r = Vec::new();
    let mut array_xmfn = Vec::new();
    let mut array_xmrn = Vec::new();
    let mut array_xmgl = Vec::new();
    let mut array_xmudnt = Vec::new();
    let mut array_rms = Vec::new();

    // Loop over radome slope values
    let mut r_val = -0.06;
    while r_val <= 0.06 + 1e-9 {
        // Initialize state variables
        let mut x1: f64 = 0.0;
        let mut x2: f64 = 0.0;
        let mut x3: f64 = 0.0;
        let mut x4: f64 = 0.0;
        let mut x5: f64 = 0.0;
        let mut x6: f64 = 0.0;
        let mut x7: f64 = 0.0;
        let mut x8: f64 = 0.0;
        let mut x9: f64 = 0.0;
        let mut x10: f64 = 0.0;
        let mut x11: f64 = 0.0;
        let mut x12: f64 = 0.0;
        let mut x13: f64 = 0.0;
        let mut x14: f64 = 0.0;
        let mut x15: f64 = 0.0;

        if miss == 1 {
            x3 = 1.0;
        } else if miss == 2 {
            x8 = 1.0;
        } else if miss == 3 {
            x6 = xnp * vc / 32.2;
        }

        let _t: f64 = 0.0;
        let mut tp: f64 = 0.00001 + tint;

        while tp <= tf - 1e-5 {
            let x1old = x1;
            let x2old = x2;
            let x3old = x3;
            let x4old = x4;
            let x5old = x5;
            let x6old = x6;
            let x7old = x7;
            let x8old = x8;
            let x9old = x9;
            let x10old = x10;
            let x11old = x11;
            let x12old = x12;
            let x13old = x13;
            let x14old = x14;
            let x15old = x15;

            // RK2 integration (two-pass)
            let mut tgo = tp + 0.00001;

            // First pass - compute derivatives
            let x1d = x2;
            let x2d = x3;
            let y1pz = (x6 / t2 + x5) / t1;
            let x3d = y1pz / (vc * tgo);
            let x4d = -y1pz;
            let x5d = -y1pz + r_val * y1pz;
            let y2pz = -xka * wi * x7;
            let x6d = -x6 / t2 + xnp * vc * xkdc * y2pz / 32.2;
            let x7d = xkr * wact * wact * x8;
            let x8d = x9 - 2.0 * zact * wact * x8;
            let y4pz = xk1 * (-32.2 * x2 - y2pz);
            let x9d = -wact * wact * x8 + x10 * waf * waf - waf * waf * y4pz / (wz * wz);
            let y3pz = xk3 * (x7d + wi * x7 + (x4 - x5) / 57.3);
            let x10d = -2.0 * zaf * waf * (x10 - y4pz / (wz * wz)) + x11 + ta * y3pz;
            let x11d = -waf * waf * (x10 - y4pz / (wz * wz)) + y4pz + y3pz;
            let x12d = x1 * x1;
            let x13d = (y1pz / (vc * tgo)).powi(2);
            let x14d = (y1pz * vc * tgo / ra).powi(2);
            let x15d = y1pz.powi(2);

            // Euler step
            x1 += h * x1d;
            x2 += h * x2d;
            x3 += h * x3d;
            x4 += h * x4d;
            x5 += h * x5d;
            x6 += h * x6d;
            x7 += h * x7d;
            x8 += h * x8d;
            x9 += h * x9d;
            x10 += h * x10d;
            x11 += h * x11d;
            x12 += h * x12d;
            x13 += h * x13d;
            x14 += h * x14d;
            x15 += h * x15d;
            tp += h;

            // Second pass - compute derivatives at new state
            tgo = tp + 0.00001;
            let x1d2 = x2;
            let x2d2 = x3;
            let y1pz2 = (x6 / t2 + x5) / t1;
            let x3d2 = y1pz2 / (vc * tgo);
            let x4d2 = -y1pz2;
            let x5d2 = -y1pz2 + r_val * y1pz2;
            let y2pz2 = -xka * wi * x7;
            let x6d2 = -x6 / t2 + xnp * vc * xkdc * y2pz2 / 32.2;
            let x7d2 = xkr * wact * wact * x8;
            let x8d2 = x9 - 2.0 * zact * wact * x8;
            let y4pz2 = xk1 * (-32.2 * x2 - y2pz2);
            let x9d2 = -wact * wact * x8 + x10 * waf * waf - waf * waf * y4pz2 / (wz * wz);
            let y3pz2 = xk3 * (x7d2 + wi * x7 + (x4 - x5) / 57.3);
            let x10d2 = -2.0 * zaf * waf * (x10 - y4pz2 / (wz * wz)) + x11 + ta * y3pz2;
            let x11d2 = -waf * waf * (x10 - y4pz2 / (wz * wz)) + y4pz2 + y3pz2;
            let x12d2 = x1 * x1;
            let x13d2 = (y1pz2 / (vc * tgo)).powi(2);
            let x14d2 = (y1pz2 * vc * tgo / ra).powi(2);
            let x15d2 = y1pz2.powi(2);

            // RK2 averaging
            x1 = 0.5 * (x1old + x1 + h * x1d2);
            x2 = 0.5 * (x2old + x2 + h * x2d2);
            x3 = 0.5 * (x3old + x3 + h * x3d2);
            x4 = 0.5 * (x4old + x4 + h * x4d2);
            x5 = 0.5 * (x5old + x5 + h * x5d2);
            x6 = 0.5 * (x6old + x6 + h * x6d2);
            x7 = 0.5 * (x7old + x7 + h * x7d2);
            x8 = 0.5 * (x8old + x8 + h * x8d2);
            x9 = 0.5 * (x9old + x9 + h * x9d2);
            x10 = 0.5 * (x10old + x10 + h * x10d2);
            x11 = 0.5 * (x11old + x11 + h * x11d2);
            x12 = 0.5 * (x12old + x12 + h * x12d2);
            x13 = 0.5 * (x13old + x13 + h * x13d2);
            x14 = 0.5 * (x14old + x14 + h * x14d2);
            x15 = 0.5 * (x15old + x15 + h * x15d2);
        }

        // Compute miss distance statistics
        let tgo = tp + 0.00001;
        let xmfn = (x15 * phifn).sqrt();
        let xmrn = (x14 * phirn).sqrt();
        let xmgl = (x13 * phigl).sqrt();
        let xmudnt = xnt * (x12 / tgo).sqrt();
        let rms_val = (xmfn.powi(2) + xmrn.powi(2) + xmgl.powi(2) + xmudnt.powi(2)).sqrt();

        array_r.push(r_val);
        array_xmfn.push(xmfn);
        array_xmrn.push(xmrn);
        array_xmgl.push(xmgl);
        array_xmudnt.push(xmudnt);
        array_rms.push(rms_val);

        r_val += 0.01;
    }

    Results {
        r: array_r,
        xmfn: array_xmfn,
        xmrn: array_xmrn,
        xmgl: array_xmgl,
        xmudnt: array_xmudnt,
        rms: array_rms,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c23l3_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.r.clone(),
        results.xmfn.clone(),
        results.xmrn.clone(),
        results.xmgl.clone(),
        results.xmudnt.clone(),
        results.rms.clone(),
    ])?;

    let plot_file = format!("{}/c23l3_miss.png", output_dir);
    let config = PlotConfig::new("Radome Slope Analysis - Miss Distance")
        .with_labels("Radome Slope", "Standard Deviation of Miss (Ft)");

    let series = vec![
        Series::new(results.r.clone(), results.xmfn.clone())
            .with_color(plotters::prelude::BLUE)
            .with_label("XMFN"),
        Series::new(results.r.clone(), results.xmrn.clone())
            .with_color(plotters::prelude::RED)
            .with_label("XMRN"),
        Series::new(results.r.clone(), results.xmgl.clone())
            .with_color(plotters::prelude::GREEN)
            .with_label("XMGL"),
        Series::new(results.r.clone(), results.xmudnt.clone())
            .with_color(plotters::prelude::MAGENTA)
            .with_label("XMUDNT"),
        Series::new(results.r.clone(), results.rms.clone())
            .with_color(plotters::prelude::BLACK)
            .with_label("RMS"),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C23L3: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c23l3_runs() {
        let results = run();
        assert_eq!(results.r.len(), 13);
    }
}
