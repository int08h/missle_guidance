//! Chapter 24, Lesson 2: Nonlinear Autopilot with Flexible Body
//!
//! Nonlinear autopilot simulation with detailed aerodynamics.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub xncg: Vec<f64>,
    pub xnlg: Vec<f64>,
    pub delnl: Vec<f64>,
    pub alf: Vec<f64>,
}

/// Run the C24L2 simulation
pub fn run() -> Results {
    // Parameters
    let xnc: f64 = 322.0;
    let xncg: f64 = xnc / 32.2;
    let alt: f64 = 0.0;
    let ts: f64 = 0.01;
    let h: f64 = 0.001;
    let vm: f64 = 3000.0;
    let wcr: f64 = 50.0;
    let zeta: f64 = 0.7;
    let tau: f64 = 0.3;
    let wact: f64 = 150.0;
    let zact: f64 = 0.7;
    let slope: f64 = 1.5;
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
    let adel = 1.33 * xlp * diam / 2.0;
    let asq = (xl - xlp) * diam;
    let xcpb = (2.0 * adel * xlp / 3.0 + asq * (xlp + 0.5 * (xl - xlp))) / (adel + asq);
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
    let _a0 = -p1 * y6;
    let alftr = (-p3 + (p3 * p3 + 4.0 * p2 * p1).sqrt()) / (2.0 * p2);
    let _deltr = -y4 * alftr / y6 - y5 * alftr * alftr / y6;
    let _alfrqd = alftr;

    let cna = 2.0 + slope * splan * alftr / sref + 8.0 * swing / (b * sref) + 8.0 * stail / (b * sref);
    let cnd = 8.0 * stail / (b * sref);
    let cmap = 2.0 * tmp4 + slope * splan * alftr * tmp3 / sref + 8.0 * swing * tmp1 / (b * sref);
    let cma = cmap + 8.0 * stail * tmp2 / (b * sref);
    let cmd = 8.0 * stail * tmp2 / (b * sref);

    let xma = q * sref * diam * cma / xiyy;
    let xmd = q * sref * diam * cmd / xiyy;
    let za = -32.2 * q * sref * cna / (wgt * vm);
    let zd = -32.2 * q * sref * cnd / (wgt * vm);

    let wz = ((xma * zd - za * xmd) / zd).sqrt();
    let b11 = za / xma;
    let b12 = -1.0 / xma;
    let xk1 = -vm * (xma * zd - xmd * za) / (1845.0 * xma);
    let _xk2 = xk1;
    let ta = xmd / (xma * zd - xmd * za);
    let xk3 = 1845.0 * xk1 / vm;

    // Autopilot gains
    let w = (tau * wcr * (1.0 + b11 / (wcr * b12)) - 1.0) / (2.0 * zeta * tau);
    let w0 = w / (tau * wcr).sqrt();
    let z0 = 0.5 * w0 * (2.0 * zeta / w + tau - 1.0 / (w0 * w0 * wcr * b12));
    let xkc = (-w0 * w0 / (wz * wz) - 1.0 + 2.0 * z0 * w0 * ta) / (1.0 - 2.0 * z0 * w0 * ta + w0 * w0 * ta * ta);
    let xka = xk3 / (xk1 * xkc);
    let xk0 = -b12 * w * w / tau;
    let _xk = xk0 / (xk1 * (1.0 + xkc));
    let wi = xkc * ta * w0 * w0 / (1.0 + xkc + w0 * w0 / (wz * wz));
    let xkr = _xk / (xka * wi);
    let xkdc = 1.0 + 1845.0 / (xka * vm);

    // State variables
    let mut t: f64 = 0.0;
    let mut s: f64 = 0.0;
    let mut thd: f64 = 0.0;
    let mut alf: f64 = 0.0;
    let mut xx: f64 = 0.0;
    let mut delnl: f64 = 0.0;
    let mut delnld: f64 = 0.0;

    let mut array_t = Vec::new();
    let mut array_xncg = Vec::new();
    let mut array_xnlg = Vec::new();
    let mut array_delnl = Vec::new();
    let mut array_alf = Vec::new();

    while t <= 0.999999 {
        s += h;

        let thdold = thd;
        let alfold = alf;
        let xxold = xx;
        let delnlold = delnl;
        let delnldold = delnld;

        // First derivative evaluation
        let cn = 2.0 * alf + 1.5 * splan * alf * alf / sref
            + 8.0 * swing * alf / (b * sref)
            + 8.0 * stail * (alf + delnl / 57.3) / (b * sref);
        let cm = 2.0 * alf * tmp4 + 1.5 * splan * alf * alf * tmp3 / sref
            + 8.0 * swing * alf * tmp1 / (b * sref)
            + 8.0 * stail * (alf + delnl / 57.3) * tmp2 / (b * sref);

        let thdd = q * sref * diam * cm / xiyy;
        let thddeg = thd * 57.3;
        let xnl = 32.2 * q * sref * cn / wgt;
        let xnlg = xnl / 32.2;
        let alfd = thd - xnl / vm;
        let xnanl = xnlg;
        let xxd = wi * (thddeg + xka * (xnanl - xncg * xkdc));
        let delcnl = xkr * (xx + thddeg);
        let delnldd = wact * wact * (delcnl - delnl - 2.0 * zact * delnld / wact);

        // Euler step
        thd += h * thdd;
        alf += h * alfd;
        xx += h * xxd;
        delnl += h * delnld;
        delnld += h * delnldd;
        t += h;

        // Second derivative for RK2
        let cn = 2.0 * alf + 1.5 * splan * alf * alf / sref
            + 8.0 * swing * alf / (b * sref)
            + 8.0 * stail * (alf + delnl / 57.3) / (b * sref);
        let cm = 2.0 * alf * tmp4 + 1.5 * splan * alf * alf * tmp3 / sref
            + 8.0 * swing * alf * tmp1 / (b * sref)
            + 8.0 * stail * (alf + delnl / 57.3) * tmp2 / (b * sref);

        let thdd = q * sref * diam * cm / xiyy;
        let thddeg = thd * 57.3;
        let xnl = 32.2 * q * sref * cn / wgt;
        let xnlg = xnl / 32.2;
        let alfd = thd - xnl / vm;
        let xnanl = xnlg;
        let xxd = wi * (thddeg + xka * (xnanl - xncg * xkdc));
        let delcnl = xkr * (xx + thddeg);
        let delnldd = wact * wact * (delcnl - delnl - 2.0 * zact * delnld / wact);

        // RK2 averaging
        thd = 0.5 * (thdold + thd + h * thdd);
        alf = 0.5 * (alfold + alf + h * alfd);
        xx = 0.5 * (xxold + xx + h * xxd);
        delnl = 0.5 * (delnlold + delnl + h * delnld);
        delnld = 0.5 * (delnldold + delnld + h * delnldd);

        if s >= ts - 0.00001 {
            s = 0.0;
            let xncg_out = xnc / 32.2;
            let xnlg_out = xnl / 32.2;
            array_t.push(t);
            array_xncg.push(xncg_out);
            array_xnlg.push(xnlg_out);
            array_delnl.push(delnl);
            array_alf.push(alf * 57.3);
        }
    }

    Results {
        time: array_t,
        xncg: array_xncg,
        xnlg: array_xnlg,
        delnl: array_delnl,
        alf: array_alf,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c24l2_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.xncg.clone(),
        results.xnlg.clone(),
        results.delnl.clone(),
        results.alf.clone(),
    ])?;

    // Acceleration plot
    let plot_file = format!("{}/c24l2_accel.png", output_dir);
    let config = PlotConfig::new("Nonlinear Autopilot - Acceleration")
        .with_labels("Time (Sec)", "Acceleration (G)");

    let series = vec![
        Series::new(results.time.clone(), results.xncg.clone())
            .with_color(plotters::prelude::RED)
            .with_label("Commanded"),
        Series::new(results.time.clone(), results.xnlg.clone())
            .with_color(plotters::prelude::BLUE)
            .with_label("Actual"),
    ];

    line_plot(&plot_file, &config, &series).ok();

    // Fin deflection plot
    let plot_file2 = format!("{}/c24l2_fin.png", output_dir);
    let config2 = PlotConfig::new("Nonlinear Autopilot - Fin Deflection")
        .with_labels("Time (Sec)", "Fin Deflection (deg)");

    let series2 = vec![
        Series::new(results.time.clone(), results.delnl.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file2, &config2, &series2).ok();

    // Angle of attack plot
    let plot_file3 = format!("{}/c24l2_aoa.png", output_dir);
    let config3 = PlotConfig::new("Nonlinear Autopilot - Angle of Attack")
        .with_labels("Time (Sec)", "Angle of Attack (deg)");

    let series3 = vec![
        Series::new(results.time.clone(), results.alf.clone())
            .with_color(plotters::prelude::GREEN),
    ];

    line_plot(&plot_file3, &config3, &series3).ok();

    println!("C24L2: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c24l2_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
