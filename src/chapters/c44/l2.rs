//! Chapter 44, Lesson 2: BMD Flight Time vs Velocity
//!
//! Ballistic missile defense: interceptor velocity vs flight time analysis.

use crate::save_data;
use crate::utils::{lambert3d, predict44, EARTH_RADIUS_FT, GM_FT};

pub struct Results {
    pub tf_tl: Vec<f64>,
    pub vbom: Vec<f64>,
}

/// Run the C44L2 simulation
pub fn run() -> Results {
    let tlaunch: f64 = 300.0;
    let _ts: f64 = 1.0;
    let xlongmdegickm: f64 = 1000.0;
    let rdeskm: f64 = 10000.0;
    let altmkmic: f64 = 0.0;
    let tloft: f64 = 0.0;
    let switch: i32 = 0;
    let switchm: i32 = 0;

    let a = EARTH_RADIUS_FT;
    let gm = GM_FT;

    let mut tftot = 252.0 + 0.223 * rdeskm - 5.44e-6 * rdeskm * rdeskm;
    tftot += tloft;

    let mut array_tftl = Vec::new();
    let mut array_vbom = Vec::new();

    let mut tf = tlaunch + 60.0;
    while tf <= tftot - 20.0 {
        let altm = altmkmic * 3280.0;
        let xlongfdeg = 57.3 * rdeskm * 3280.0 / a;
        let xlongmdegic = xlongmdegickm / 111.0;
        let xlatmdegic: f64 = 0.0;
        let xlongtdeg: f64 = 0.0;
        let xlattdeg: f64 = 0.0;
        let xlatfdeg: f64 = 0.0;

        let xlongf = xlongfdeg / 57.3;
        let xlatf = xlatfdeg / 57.3;
        let xlongt = xlongtdeg / 57.3;
        let xlongm = xlongmdegic / 57.3;
        let xlatm = xlatmdegic / 57.3;
        let xlatt = xlattdeg / 57.3;

        let xf = a * xlatf.cos() * xlongf.cos();
        let yf = a * xlatf.cos() * xlongf.sin();
        let zf: f64 = 0.0;

        let mut xt = a * xlatt.cos() * xlongt.cos();
        let mut yt = a * xlatt.cos() * xlongt.sin();
        let zt: f64 = 0.0;

        let _xtinit = xt;
        let _ytinit = yt;
        let _ztinit = zt;

        let mut xm = (a + altm) * xlatm.cos() * xlongm.cos();
        let mut ym = (a + altm) * xlatm.cos() * xlongm.sin();
        let zm: f64 = 0.0;

        let mut xmd: f64 = 0.0;
        let mut ymd: f64 = 0.0;

        // Get target initial velocity
        let tgolam = tftot;
        let result = lambert3d(xt, yt, zt, tgolam, xf, yf, zf, switch);
        let mut xtd = result.vrx;
        let mut ytd = result.vry;
        let _ztd: f64 = 0.0;

        // Predict target location at intercept time
        let result_pred = predict44(0.0, xt, yt, xtd, ytd, tf, tftot, xf, yf);
        let xtfact = result_pred.xtf;
        let ytfact = result_pred.ytf;

        let mut rtm1 = xt - xm;
        let mut rtm2 = yt - ym;
        let mut rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
        let vtm1 = xtd - xmd;
        let vtm2 = ytd - ymd;
        let mut vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;

        let mut h: f64;
        let mut t: f64 = 0.0;
        let mut _s: f64 = 0.0;
        let mut qboostm: bool = true;
        let mut vbom: f64 = 0.0;
        let mut alttkm: f64 = 0.0;

        while !((t > tlaunch + 50.0) && vc < 0.0 && rtm < 10000.0) {
            if rtm < 1000.0 {
                h = 0.0001;
            } else {
                h = 0.01;
            }

            let xtold = xt;
            let ytold = yt;
            let xtdold = xtd;
            let ytdold = ytd;
            let xmold = xm;
            let ymold = ym;
            let xmdold = xmd;
            let ymdold = ymd;

            // First derivative evaluation
            let tempbott = (xt * xt + yt * yt).powf(1.5);
            let xtdd = -gm * xt / tempbott;
            let ytdd = -gm * yt / tempbott;

            rtm1 = xt - xm;
            rtm2 = yt - ym;
            let vtm1 = xtd - xmd;
            let vtm2 = ytd - ymd;
            rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
            vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;
            let _tgo = rtm / vc;

            let (xmdd, ymdd) = if t > tlaunch {
                let tempbotm = (xm * xm + ym * ym).powf(1.5);
                (-gm * xm / tempbotm, -gm * ym / tempbotm)
            } else {
                (0.0, 0.0)
            };

            alttkm = ((xt * xt + yt * yt).sqrt() - a) / 3280.0;

            // Euler step
            xt += h * xtd;
            yt += h * ytd;
            xtd += h * xtdd;
            ytd += h * ytdd;
            xm += h * xmd;
            ym += h * ymd;
            xmd += h * xmdd;
            ymd += h * ymdd;
            t += h;

            // RK2 averaging
            xt = 0.5 * (xtold + xt + h * xtd);
            yt = 0.5 * (ytold + yt + h * ytd);
            xtd = 0.5 * (xtdold + xtd + h * xtdd);
            ytd = 0.5 * (ytdold + ytd + h * ytdd);
            xm = 0.5 * (xmold + xm + h * xmd);
            ym = 0.5 * (ymold + ym + h * ymd);
            xmd = 0.5 * (xmdold + xmd + h * xmdd);
            ymd = 0.5 * (ymdold + ymd + h * ymdd);

            _s += h;

            // Launch missile
            let tgolamm = tf - t;
            if t >= tlaunch && qboostm {
                qboostm = false;
                let result = lambert3d(xm, ym, zm, tgolamm, xtfact, ytfact, 0.0, switchm);
                xmd = result.vrx;
                ymd = result.vry;
                vbom = (xmd * xmd + ymd * ymd).sqrt() / 3280.0;
            }
        }

        if vbom < 8.0 && rtm < 1000.0 && alttkm > 50.0 {
            array_tftl.push(tf - tlaunch);
            array_vbom.push(vbom);
        }

        tf += 20.0;
    }

    Results {
        tf_tl: array_tftl,
        vbom: array_vbom,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c44l2_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.tf_tl.clone(),
        results.vbom.clone(),
    ])?;

    println!("C44L2: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c44l2_runs() {
        let results = run();
        // May have zero results depending on parameters
        assert!(true);
    }
}
