//! Chapter 21, Lesson 2: Missile Aerodynamics with Transfer Function
//!
//! Computes aerodynamic response with linearized transfer function model,
//! including acceleration and pitch rate outputs.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub xnl: Vec<f64>,   // Missile Acceleration (G)
    pub thd: Vec<f64>,   // Pitch rate
}

/// Run the C21L2 simulation
pub fn run() -> Results {
    let vm: f64 = 3000.0;
    let del: f64 = 5.0 / 57.3;  // Fin deflection (rad)
    let alt: f64 = 0.0;
    let a: f64 = 1000.0;        // Speed of sound
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

    // Trim condition calculation
    let y1 = 2.0 * tmp4 + 8.0 * swing * tmp1 / (b * sref) + 8.0 * stail * tmp2 / (b * sref);
    let y2 = 1.5 * splan * tmp3 / sref;
    let y3 = 8.0 * stail * tmp2 * del / (b * sref);
    let alftr = (-y1 - (y1 * y1 - 4.0 * y2 * y3).sqrt()) / (2.0 * y2);

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
    let ta = xmd / (xma * zd - xmd * za);
    let xk3 = 1845.0 * xk1 / vm;

    let mut e: f64 = 0.0;
    let mut ed: f64 = 0.0;
    let mut t: f64 = 0.0;
    let h: f64 = 0.0025;
    let mut s: f64 = 0.0;

    let mut array_t = Vec::new();
    let mut array_xnl = Vec::new();
    let mut array_thd = Vec::new();

    let mut xnl: f64;
    let mut thd: f64;

    while t < 1.99999 {
        let eold = e;
        let edold = ed;

        // First derivative evaluation
        let edd = waf * waf * (del * 57.3 - e - 2.0 * zaf * ed / waf);
        let _xnl = xk1 * (e - edd / (wz * wz));
        let _thd = xk3 * (e + ta * ed);

        // Euler step
        e += h * ed;
        ed += h * edd;
        t += h;

        // Second derivative for RK2
        let edd = waf * waf * (del * 57.3 - e - 2.0 * zaf * ed / waf);
        xnl = xk1 * (e - edd / (wz * wz));
        thd = xk3 * (e + ta * ed);

        // RK2 averaging
        e = 0.5 * (eold + e + h * ed);
        ed = 0.5 * (edold + ed + h * edd);

        s += h;
        if s >= 0.0099999 {
            s = 0.0;
            array_t.push(t);
            array_xnl.push(xnl);
            array_thd.push(thd);
        }
    }

    Results {
        time: array_t,
        xnl: array_xnl,
        thd: array_thd,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c21l2_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.xnl.clone(),
        results.thd.clone(),
    ])?;

    let plot_file = format!("{}/c21l2_accel.png", output_dir);
    let config = PlotConfig::new("Missile Aerodynamics - Transfer Function Response")
        .with_labels("Time (Sec)", "Missile Acceleration (G)");

    let series = vec![
        Series::new(results.time.clone(), results.xnl.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C21L2: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c21l2_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
