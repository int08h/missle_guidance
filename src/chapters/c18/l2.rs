//! Chapter 18, Lesson 2: Strategic Intercept with Boosting Target
//!
//! PN guidance for intercepting a boosting two-stage rocket target.
//!
//! MATLAB BUG: Calls `predictb` and `lambertpz`. The lambertpz function file
//! is `lamberpz.m` and the function is defined as `lambert` not `lambertpz`.
//! Output: 8 columns [T, DISTNMT, ALTNMT, DISTNMM, ALTNMM, XNCG, DELV, ATPLOSG]

use crate::save_data;
use crate::utils::lambert2d::lambert2d;
use crate::utils::predict::predict_b;
use crate::utils::constants::G_ACCEL;

pub struct Results {
    pub time: Vec<f64>,
    pub distnmt: Vec<f64>,
    pub altnmt: Vec<f64>,
    pub distnmm: Vec<f64>,
    pub altnmm: Vec<f64>,
    pub xncg: Vec<f64>,
    pub delv: Vec<f64>,
    pub atplosg: Vec<f64>,
}

/// Run the C18L2 simulation
pub fn run() -> Results {
    let left: i32 = 0;
    let xnclim: f64 = 644.0;
    let a: f64 = 2.0926e7;
    let gm: f64 = 1.4077e16;
    let degrad = 360.0 / (2.0 * std::f64::consts::PI);
    let xnp: f64 = 3.0;
    let prederr: f64 = 0.0;

    // Two-stage rocket parameters
    let xisp1: f64 = 250.0;
    let xisp2: f64 = 250.0;
    let xmf1: f64 = 0.85;
    let xmf2: f64 = 0.85;
    let wpay: f64 = 100.0;
    let delv_total: f64 = 20000.0;
    let delv1 = 0.3333 * delv_total;
    let delv2 = 0.6667 * delv_total;
    let amax1: f64 = 20.0;
    let amax2: f64 = 20.0;
    let xkickdeg: f64 = 80.0;

    // Stage 2 calculations
    let top2 = wpay * ((delv2 / (xisp2 * G_ACCEL)).exp() - 1.0);
    let bot2 = 1.0 / xmf2 - ((1.0 - xmf2) / xmf2) * (delv2 / (xisp2 * G_ACCEL)).exp();
    let wp2 = top2 / bot2;
    let ws2 = wp2 * (1.0 - xmf2) / xmf2;
    let wtot2 = wp2 + ws2 + wpay;
    let trst2 = amax2 * (wpay + ws2);
    let tb2 = xisp2 * wp2 / trst2;

    // Stage 1 calculations
    let top1 = wtot2 * ((delv1 / (xisp1 * G_ACCEL)).exp() - 1.0);
    let bot1 = 1.0 / xmf1 - ((1.0 - xmf1) / xmf1) * (delv1 / (xisp1 * G_ACCEL)).exp();
    let wp1 = top1 / bot1;
    let ws1 = wp1 * (1.0 - xmf1) / xmf1;
    let wtot = wp1 + ws1 + wtot2;
    let trst1 = amax1 * (wtot2 + ws1);
    let tb1 = xisp1 * wp1 / trst1;

    let xlongmdeg: f64 = 85.0;
    let xlongtdeg: f64 = 90.0;
    let altnmtic: f64 = 0.0;
    let altnmmic: f64 = 0.0;
    let tf: f64 = 50.0;

    let altt_init = altnmtic * 6076.0;
    let altm_init = altnmmic * 6076.0;
    let mut s: f64 = 0.0;
    let xlongm = xlongmdeg / degrad;
    let xlongt = xlongtdeg / degrad;

    let mut xm = (a + altm_init) * xlongm.cos();
    let mut ym = (a + altm_init) * xlongm.sin();
    let mut xt = (a + altt_init) * xlongt.cos();
    let mut yt = (a + altt_init) * xlongt.sin();
    let xfirstt = xt;
    let yfirstt = yt;

    let (x1t_init, y1t_init) = if left == 1 {
        (
            (std::f64::consts::PI / 2.0 - xkickdeg / degrad + xlongt).cos(),
            (std::f64::consts::PI / 2.0 - xkickdeg / degrad + xlongt).sin(),
        )
    } else {
        (
            (-std::f64::consts::PI / 2.0 + xkickdeg / degrad + xlongt).cos(),
            (-std::f64::consts::PI / 2.0 + xkickdeg / degrad + xlongt).sin(),
        )
    };

    let mut t: f64 = 0.0;

    // Predict target position at TF using two-stage boost
    let pred = predict_b(tf, xt, yt, x1t_init, y1t_init, wp1, wtot, tb1, trst1, tb2, wp2, wtot2, trst2, wpay);
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
    let mut array_xncg = Vec::new();
    let mut array_delv = Vec::new();
    let mut array_atplosg = Vec::new();

    let mut xnc: f64;
    let mut atplos: f64;

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
        let (wgt, trst) = if t < tb1 {
            (-wp1 * t / tb1 + wtot, trst1)
        } else if t < tb1 + tb2 {
            (-wp2 * t / tb2 + wtot2 + wp2 * tb1 / tb2, trst2)
        } else {
            (wpay, 0.0)
        };

        let at = G_ACCEL * trst / wgt;
        let vel = (x1t * x1t + y1t * y1t).sqrt();
        let axt = at * x1t / vel;
        let ayt = at * y1t / vel;

        let tembott = (xt * xt + yt * yt).powf(1.5);
        let mut x1dt = -gm * xt / tembott + axt;
        let mut y1dt = -gm * yt / tembott + ayt;
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
        // atplos computed in second evaluation for RK2

        xnc = if t > 25.0 { xnp * vc * xlamd } else { 0.0 };
        xnc = xnc.clamp(-xnclim, xnclim);

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
        let (wgt, trst) = if t < tb1 {
            (-wp1 * t / tb1 + wtot, trst1)
        } else if t < tb1 + tb2 {
            (-wp2 * t / tb2 + wtot2 + wp2 * tb1 / tb2, trst2)
        } else {
            (wpay, 0.0)
        };

        let at = G_ACCEL * trst / wgt;
        let vel = (x1t * x1t + y1t * y1t).sqrt();
        let axt = at * x1t / vel;
        let ayt = at * y1t / vel;

        let tembott = (xt * xt + yt * yt).powf(1.5);
        x1dt = -gm * xt / tembott + axt;
        y1dt = -gm * yt / tembott + ayt;
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
        atplos = y1dt * xlam.cos() - x1dt * xlam.sin();

        xnc = if t > 25.0 { xnp * vc * xlamd } else { 0.0 };
        xnc = xnc.clamp(-xnclim, xnclim);

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
            let xncg = xnc / G_ACCEL;
            let r = (xm * xm + ym * ym).sqrt();
            let cbeta = (xm * xfirstt + ym * yfirstt) / (r * rf);
            let beta = cbeta.acos();
            let distnmm = a * beta / 6076.0;
            let atplosg = atplos / G_ACCEL;

            array_t.push(t);
            array_distnmt.push(distnmt);
            array_altnmt.push(altnmt);
            array_distnmm.push(distnmm);
            array_altnmm.push(altnmm);
            array_xncg.push(xncg);
            array_delv.push(delv);
            array_atplosg.push(atplosg);
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
        xncg: array_xncg,
        delv: array_delv,
        atplosg: array_atplosg,
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
        results.xncg.clone(),
        results.delv.clone(),
        results.atplosg.clone(),
    ];
    let data_file = format!("{}/c18l2_datfil.txt", output_dir);
    save_data(&data_file, &data)?;

    println!("C18L2: Simulation complete");
    println!("  Output: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c18l2_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
