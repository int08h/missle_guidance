//! Chapter 18, Lesson 1: Strategic Intercept with Ballistic Target
//!
//! PN guidance for intercepting a ballistic target on a spherical Earth.
//!
//! MATLAB BUG: Calls `predictpz` and `lambertpz`. The lambertpz function file
//! is `lamberpz.m` and the function is defined as `lambert` not `lambertpz`.
//! Output: 5 columns [T, DISTNMT, ALTNMT, DISTNMM, ALTNMM]

use crate::save_data;
use crate::utils::lambert2d::lambert2d;
use crate::utils::predict::predict_pz;

pub struct Results {
    pub time: Vec<f64>,
    pub distnmt: Vec<f64>,
    pub altnmt: Vec<f64>,
    pub distnmm: Vec<f64>,
    pub altnmm: Vec<f64>,
}

/// Run the C18L1 simulation
pub fn run() -> Results {
    let xlongmdeg: f64 = 45.0;
    let xlongtdeg: f64 = 90.0;
    let altnmtic: f64 = 0.0;
    let altnmmic: f64 = 0.0;
    let tf: f64 = 500.0;
    let gamdegt: f64 = 23.0;
    let a: f64 = 2.0926e7;
    let gm: f64 = 1.4077e16;
    let degrad = 360.0 / (2.0 * std::f64::consts::PI);
    let xnp: f64 = 3.0;
    let prederr: f64 = 0.0;

    let gamt = gamdegt / 57.3;
    let distnmt_ic: f64 = 6000.0;
    let phit = distnmt_ic * 6076.0 / a;
    let altt_init = altnmtic * 6076.0;
    let altm_init = altnmmic * 6076.0;
    let r0t = a + altt_init;
    let top = gm * (1.0 - phit.cos());
    let temp = r0t * gamt.cos() / a - (phit + gamt).cos();
    let bot = r0t * gamt.cos() * temp;
    let vt = (top / bot).sqrt();

    let xlongm = xlongmdeg / degrad;
    let xlongt = xlongtdeg / degrad;

    let (x1t_init, y1t_init) = if xlongm > xlongt {
        (
            vt * (std::f64::consts::PI / 2.0 - gamt + xlongt).cos(),
            vt * (std::f64::consts::PI / 2.0 - gamt + xlongt).sin(),
        )
    } else {
        (
            vt * (-std::f64::consts::PI / 2.0 + gamt + xlongt).cos(),
            vt * (-std::f64::consts::PI / 2.0 + gamt + xlongt).sin(),
        )
    };

    let mut s: f64 = 0.0;
    let mut xm = (a + altm_init) * xlongm.cos();
    let mut ym = (a + altm_init) * xlongm.sin();
    let mut xt = (a + altt_init) * xlongt.cos();
    let mut yt = (a + altt_init) * xlongt.sin();
    let xfirstt = xt;
    let yfirstt = yt;
    let mut t: f64 = 0.0;

    // Predict target position at TF
    let pred = predict_pz(tf, xt, yt, x1t_init, y1t_init);
    let xtf = pred.xtf;
    let ytf = pred.ytf + prederr;

    // Lambert solver for initial missile velocity
    let lamb = lambert2d(xm, ym, tf, xtf, ytf, xlongm, xlongt);
    let mut x1m = lamb.vrx;
    let mut y1m = lamb.vry;
    let mut x1t = x1t_init;
    let mut y1t = y1t_init;

    let mut rtm1 = xt - xm;
    let mut rtm2 = yt - ym;
    let mut rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
    let vtm1 = x1t - x1m;
    let vtm2 = y1t - y1m;
    let mut vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;
    let mut delv: f64 = 0.0;

    let mut array_t = Vec::new();
    let mut array_distnmt = Vec::new();
    let mut array_altnmt = Vec::new();
    let mut array_distnmm = Vec::new();
    let mut array_altnmm = Vec::new();

    while vc >= 0.0 {
        let tgo = rtm / vc;
        let h = if tgo > 0.1 { 0.01 } else { 0.0001 };

        let xoldt = xt;
        let yoldt = yt;
        let x1oldt = x1t;
        let y1oldt = y1t;
        let xoldm = xm;
        let yoldm = ym;
        let x1oldm = x1m;
        let y1oldm = y1m;
        let delvold = delv;

        // First Euler step
        let tembott = (xt * xt + yt * yt).powf(1.5);
        let mut x1dt = -gm * xt / tembott;
        let mut y1dt = -gm * yt / tembott;
        let mut xdt = x1t;
        let mut ydt = y1t;

        rtm1 = xt - xm;
        rtm2 = yt - ym;
        rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
        let vtm1 = x1t - x1m;
        let vtm2 = y1t - y1m;
        vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;

        let xlam = rtm2.atan2(rtm1);
        let xlamd = (rtm1 * vtm2 - rtm2 * vtm1) / (rtm * rtm);
        let xnc = xnp * vc * xlamd;
        let mut delvd = xnc.abs();
        let am1 = -xnc * xlam.sin();
        let am2 = xnc * xlam.cos();
        let tembotm = (xm * xm + ym * ym).powf(1.5);
        let mut x1dm = -gm * xm / tembotm + am1;
        let mut y1dm = -gm * ym / tembotm + am2;
        let mut xdm = x1m;
        let mut ydm = y1m;

        // Euler step
        xt += h * xdt;
        yt += h * ydt;
        x1t += h * x1dt;
        y1t += h * y1dt;
        xm += h * xdm;
        ym += h * ydm;
        x1m += h * x1dm;
        y1m += h * y1dm;
        delv += h * delvd;
        t += h;

        // Second evaluation for RK2
        let tembott = (xt * xt + yt * yt).powf(1.5);
        x1dt = -gm * xt / tembott;
        y1dt = -gm * yt / tembott;
        xdt = x1t;
        ydt = y1t;

        rtm1 = xt - xm;
        rtm2 = yt - ym;
        rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
        let vtm1 = x1t - x1m;
        let vtm2 = y1t - y1m;
        vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;

        let xlam = rtm2.atan2(rtm1);
        let xlamd = (rtm1 * vtm2 - rtm2 * vtm1) / (rtm * rtm);
        let xnc = xnp * vc * xlamd;
        delvd = xnc.abs();
        let am1 = -xnc * xlam.sin();
        let am2 = xnc * xlam.cos();
        let tembotm = (xm * xm + ym * ym).powf(1.5);
        x1dm = -gm * xm / tembotm + am1;
        y1dm = -gm * ym / tembotm + am2;
        xdm = x1m;
        ydm = y1m;

        // RK2 averaging
        xt = (xoldt + xt) / 2.0 + 0.5 * h * xdt;
        yt = (yoldt + yt) / 2.0 + 0.5 * h * ydt;
        x1t = (x1oldt + x1t) / 2.0 + 0.5 * h * x1dt;
        y1t = (y1oldt + y1t) / 2.0 + 0.5 * h * y1dt;
        xm = (xoldm + xm) / 2.0 + 0.5 * h * xdm;
        ym = (yoldm + ym) / 2.0 + 0.5 * h * ydm;
        x1m = (x1oldm + x1m) / 2.0 + 0.5 * h * x1dm;
        y1m = (y1oldm + y1m) / 2.0 + 0.5 * h * y1dm;
        delv = (delvold + delv) / 2.0 + 0.5 * h * delvd;

        let altt = (xt * xt + yt * yt).sqrt() - a;
        let altm = (xm * xm + ym * ym).sqrt() - a;

        s += h;
        if s >= 0.99999 {
            s = 0.0;
            let altnmt = altt / 6076.0;
            let r = (xt * xt + yt * yt).sqrt();
            let rf = (xfirstt * xfirstt + yfirstt * yfirstt).sqrt();
            let cbeta = (xt * xfirstt + yt * yfirstt) / (r * rf);
            let beta = cbeta.acos();
            let distnmt = a * beta / 6076.0;
            let altnmm = altm / 6076.0;
            let r = (xm * xm + ym * ym).sqrt();
            let cbeta = (xm * xfirstt + ym * yfirstt) / (r * rf);
            let beta = cbeta.acos();
            let distnmm = a * beta / 6076.0;

            array_t.push(t);
            array_distnmt.push(distnmt);
            array_altnmt.push(altnmt);
            array_distnmm.push(distnmm);
            array_altnmm.push(altnmm);
        }
    }

    println!("RTM = {:.6e}", rtm);
    println!("DELV = {:.6e}", delv);

    Results {
        time: array_t,
        distnmt: array_distnmt,
        altnmt: array_altnmt,
        distnmm: array_distnmm,
        altnmm: array_altnmm,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data: Vec<Vec<f64>> = vec![
        results.time.clone(),
        results.distnmt.clone(),
        results.altnmt.clone(),
        results.distnmm.clone(),
        results.altnmm.clone(),
    ];
    let data_file = format!("{}/c18l1_datfil.txt", output_dir);
    save_data(&data_file, &data)?;

    println!("C18L1: Simulation complete");
    println!("  Output: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c18l1_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
