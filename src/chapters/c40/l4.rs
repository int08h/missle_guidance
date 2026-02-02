//! Chapter 40, Lesson 4: Strategic Intercept Simulation
//!
//! Ballistic missile intercept with predictive guidance.

use crate::save_data;
use crate::utils::{lambert3d, kepler1, StateVector, distance3dkm, EARTH_RADIUS_FT, GM_FT};

pub struct Results {
    pub time: Vec<f64>,
    pub dist_rt_nm: Vec<f64>,
    pub alt_nm: Vec<f64>,
    pub dist_rm_nm: Vec<f64>,
    pub alt_m_nm: Vec<f64>,
    pub axm_guid_g: Vec<f64>,
}

/// Run the C40L4 simulation
pub fn run() -> Results {
    let switch1 = 0;
    let switchm = 0;
    let tlaunch: f64 = 200.0;
    let xlongtdeg: f64 = 7.42;
    let xlattdeg: f64 = 43.75;
    let xlatfdeg: f64 = 36.175;
    let xlongfdeg: f64 = -115.136;
    let xlongmdeg_ic: f64 = -74.423;
    let xlatmdeg_ic: f64 = 39.364;
    let prederr: f64 = 10.0 * 6076.0;
    let xnclim: f64 = 161.0;
    let tftot: f64 = 2000.0;
    let tf: f64 = 1000.0;
    let a = EARTH_RADIUS_FT;
    let gm = GM_FT;
    let w: f64 = -6.283185 / 86400.0;
    let qboostm_init = true;

    let xlongf = xlongfdeg / 57.3 - w * tf;
    let xlatf = xlatfdeg / 57.3;
    let xlongt = xlongtdeg / 57.3;
    let xlatt = xlattdeg / 57.3;
    let xlongm = xlongmdeg_ic / 57.3;
    let xlatm = xlatmdeg_ic / 57.3;

    let xf = a * xlatf.cos() * xlongf.cos();
    let yf = a * xlatf.cos() * xlongf.sin();
    let zf = a * xlatf.sin();

    let mut xt = a * xlatt.cos() * xlongt.cos();
    let mut yt = a * xlatt.cos() * xlongt.sin();
    let mut zt = a * xlatt.sin();

    let xtinit = xt;
    let ytinit = yt;
    let ztinit = zt;

    // Get target initial velocity from Lambert
    let result = lambert3d(xt, yt, zt, tftot, xf, yf, zf, switch1);
    let mut xtd = result.vrx;
    let mut ytd = result.vry;
    let mut ztd = result.vrz;

    let mut xm = a * xlatm.cos() * xlongm.cos();
    let mut ym = a * xlatm.cos() * xlongm.sin();
    let mut zm = a * xlatm.sin();
    let mut xmd = a * w * xlatm.cos() * xlongm.sin();
    let mut ymd = -a * w * xlatm.cos() * xlongm.cos();
    let mut zmd: f64 = 0.0;

    let mut rtm1 = xt - xm;
    let mut rtm2 = yt - ym;
    let mut rtm3 = zt - zm;
    let mut rtm = (rtm1 * rtm1 + rtm2 * rtm2 + rtm3 * rtm3).sqrt();

    let mut vtm1 = xtd - xmd;
    let mut vtm2 = ytd - ymd;
    let mut vtm3 = ztd - zmd;
    let mut vc = -(rtm1 * vtm1 + rtm2 * vtm2 + rtm3 * vtm3) / rtm;

    // Predict target position at intercept time using Kepler
    let t0: f64 = 0.0;
    let t1 = tf;
    let x0 = StateVector::new(xt / 3280.0, yt / 3280.0, zt / 3280.0,
                               xtd / 3280.0, ytd / 3280.0, ztd / 3280.0);
    let x1 = kepler1(&x0, t0, t1);
    let xtf = x1.x * 3280.0 + prederr;
    let ytf = x1.y * 3280.0;
    let ztf = x1.z * 3280.0;

    let mut t = 0.0;
    let mut s = 0.0;
    let mut h: f64;
    let mut delv: f64 = 0.0;
    let mut axmguid: f64 = 0.0;
    let mut aymguid: f64 = 0.0;
    let mut azmguid: f64 = 0.0;
    let mut qboostm = qboostm_init;

    let mut array_t = Vec::new();
    let mut array_distrtnm = Vec::new();
    let mut array_altnm = Vec::new();
    let mut array_distrmnm = Vec::new();
    let mut array_altmnm = Vec::new();
    let mut array_axmguidg = Vec::new();

    let _count = 0;

    while vc > 0.0 {
        if rtm > 5000.0 {
            h = 0.01;
        } else {
            h = 0.00001;
        }

        let xtold = xt;
        let ytold = yt;
        let ztold = zt;
        let xtdold = xtd;
        let ytdold = ytd;
        let ztdold = ztd;
        let xmold = xm;
        let ymold = ym;
        let zmold = zm;
        let xmdold = xmd;
        let ymdold = ymd;
        let zmdold = zmd;
        let delvold = delv;

        // First derivative evaluation (pass 1)
        let mut tempbott = (xt * xt + yt * yt + zt * zt).powf(1.5);
        let mut xtdd = -gm * xt / tempbott;
        let mut ytdd = -gm * yt / tempbott;
        let mut ztdd = -gm * zt / tempbott;

        rtm1 = xt - xm;
        rtm2 = yt - ym;
        rtm3 = zt - zm;
        vtm1 = xtd - xmd;
        vtm2 = ytd - ymd;
        vtm3 = ztd - zmd;
        rtm = (rtm1 * rtm1 + rtm2 * rtm2 + rtm3 * rtm3).sqrt();
        vc = -(rtm1 * vtm1 + rtm2 * vtm2 + rtm3 * vtm3) / rtm;
        let mut tgo = rtm / vc;

        // ZEM guidance (pass 1)
        if tgo < 200.0 && t > (tlaunch + 50.0) {
            let zem1 = rtm1 + vtm1 * tgo;
            let zem2 = rtm2 + vtm2 * tgo;
            let zem3 = rtm3 + vtm3 * tgo;
            let zemdotrtm = (zem1 * rtm1 + zem2 * rtm2 + zem3 * rtm3) / rtm;
            let zemper1 = zem1 - zemdotrtm * rtm1 / rtm;
            let zemper2 = zem2 - zemdotrtm * rtm2 / rtm;
            let zemper3 = zem3 - zemdotrtm * rtm3 / rtm;
            axmguid = 3.0 * zemper1 / (tgo * tgo);
            aymguid = 3.0 * zemper2 / (tgo * tgo);
            azmguid = 3.0 * zemper3 / (tgo * tgo);
        } else {
            axmguid = 0.0;
            aymguid = 0.0;
            azmguid = 0.0;
        }

        // Clamp guidance commands
        if axmguid > xnclim {
            axmguid = xnclim;
        } else if axmguid < -xnclim {
            axmguid = -xnclim;
        }
        if aymguid > xnclim {
            aymguid = xnclim;
        } else if aymguid < -xnclim {
            aymguid = -xnclim;
        }
        if azmguid > xnclim {
            azmguid = xnclim;
        } else if azmguid < -xnclim {
            azmguid = -xnclim;
        }

        let (mut xmdd, mut ymdd, mut zmdd) = if t > tlaunch {
            let tempbotm = (xm * xm + ym * ym + zm * zm).powf(1.5);
            (-gm * xm / tempbotm + axmguid,
             -gm * ym / tempbotm + aymguid,
             -gm * zm / tempbotm + azmguid)
        } else {
            (0.0, 0.0, 0.0)
        };

        let mut delvd = (axmguid * axmguid + aymguid * aymguid + azmguid * azmguid).sqrt();

        // Euler step
        xt += h * xtd;
        yt += h * ytd;
        zt += h * ztd;
        xtd += h * xtdd;
        ytd += h * ytdd;
        ztd += h * ztdd;

        if t < tlaunch && qboostm {
            // Missile stays on rotating Earth before launch
        } else {
            xm += h * xmd;
            ym += h * ymd;
            zm += h * zmd;
            xmd += h * xmdd;
            ymd += h * ymdd;
            zmd += h * zmdd;
        }

        delv += h * delvd;
        t += h;

        // Second derivative evaluation (pass 2) - at Euler-stepped state
        tempbott = (xt * xt + yt * yt + zt * zt).powf(1.5);
        xtdd = -gm * xt / tempbott;
        ytdd = -gm * yt / tempbott;
        ztdd = -gm * zt / tempbott;
        let _altnm = ((xt * xt + yt * yt + zt * zt).sqrt() - a) / 3280.0;

        rtm1 = xt - xm;
        rtm2 = yt - ym;
        rtm3 = zt - zm;
        vtm1 = xtd - xmd;
        vtm2 = ytd - ymd;
        vtm3 = ztd - zmd;
        rtm = (rtm1 * rtm1 + rtm2 * rtm2 + rtm3 * rtm3).sqrt();
        vc = -(rtm1 * vtm1 + rtm2 * vtm2 + rtm3 * vtm3) / rtm;
        tgo = rtm / vc;

        // ZEM guidance (pass 2)
        if tgo < 200.0 && t > (tlaunch + 50.0) {
            let zem1 = rtm1 + vtm1 * tgo;
            let zem2 = rtm2 + vtm2 * tgo;
            let zem3 = rtm3 + vtm3 * tgo;
            let zemdotrtm = (zem1 * rtm1 + zem2 * rtm2 + zem3 * rtm3) / rtm;
            let zemper1 = zem1 - zemdotrtm * rtm1 / rtm;
            let zemper2 = zem2 - zemdotrtm * rtm2 / rtm;
            let zemper3 = zem3 - zemdotrtm * rtm3 / rtm;
            axmguid = 3.0 * zemper1 / (tgo * tgo);
            aymguid = 3.0 * zemper2 / (tgo * tgo);
            azmguid = 3.0 * zemper3 / (tgo * tgo);
        } else {
            axmguid = 0.0;
            aymguid = 0.0;
            azmguid = 0.0;
        }

        // Clamp guidance commands (pass 2)
        if axmguid > xnclim {
            axmguid = xnclim;
        } else if axmguid < -xnclim {
            axmguid = -xnclim;
        }
        if aymguid > xnclim {
            aymguid = xnclim;
        } else if aymguid < -xnclim {
            aymguid = -xnclim;
        }
        if azmguid > xnclim {
            azmguid = xnclim;
        } else if azmguid < -xnclim {
            azmguid = -xnclim;
        }

        if t > tlaunch {
            let tempbotm = (xm * xm + ym * ym + zm * zm).powf(1.5);
            xmdd = -gm * xm / tempbotm + axmguid;
            ymdd = -gm * ym / tempbotm + aymguid;
            zmdd = -gm * zm / tempbotm + azmguid;
        } else {
            xmdd = 0.0;
            ymdd = 0.0;
            zmdd = 0.0;
        }

        delvd = (axmguid * axmguid + aymguid * aymguid + azmguid * azmguid).sqrt();

        // RK2 averaging (using pass 2 derivatives)
        xt = 0.5 * (xtold + xt + h * xtd);
        yt = 0.5 * (ytold + yt + h * ytd);
        zt = 0.5 * (ztold + zt + h * ztd);
        xtd = 0.5 * (xtdold + xtd + h * xtdd);
        ytd = 0.5 * (ytdold + ytd + h * ytdd);
        ztd = 0.5 * (ztdold + ztd + h * ztdd);

        if t < tlaunch && qboostm {
            xm = a * (xlatmdeg_ic / 57.3).cos() * ((xlongmdeg_ic / 57.3) - w * t).cos();
            ym = a * (xlatmdeg_ic / 57.3).cos() * ((xlongmdeg_ic / 57.3) - w * t).sin();
            zm = a * (xlatmdeg_ic / 57.3).sin();
            xmd = a * w * (xlatmdeg_ic / 57.3).cos() * ((xlongmdeg_ic / 57.3) - w * t).sin();
            ymd = -a * w * (xlatmdeg_ic / 57.3).cos() * ((xlongmdeg_ic / 57.3) - w * t).cos();
            zmd = 0.0;
        } else {
            xm = 0.5 * (xmold + xm + h * xmd);
            ym = 0.5 * (ymold + ym + h * ymd);
            zm = 0.5 * (zmold + zm + h * zmd);
            xmd = 0.5 * (xmdold + xmd + h * xmdd);
            ymd = 0.5 * (ymdold + ymd + h * ymdd);
            zmd = 0.5 * (zmdold + zmd + h * zmdd);
        }
        delv = 0.5 * (delvold + delv + h * delvd);

        // Launch missile
        let tgom = tf - t;
        if t >= tlaunch && qboostm {
            let result = lambert3d(xm, ym, zm, tgom, xtf, ytf, ztf, switchm);
            xmd = result.vrx;
            ymd = result.vry;
            zmd = result.vrz;
            qboostm = false;
        }

        s += h;
        if s >= 9.9999 {
            s = 0.0;
            let xte = xt * (w * t).cos() - yt * (w * t).sin();
            let yte = xt * (w * t).sin() + yt * (w * t).cos();
            let zte = zt;

            let distrtnm = distance3dkm(xte, yte, zte, xtinit, ytinit, ztinit);
            let altnm = ((xt * xt + yt * yt + zt * zt).sqrt() - a) / 3280.0;

            let xme = xm * (w * t).cos() - ym * (w * t).sin();
            let yme = xm * (w * t).sin() + ym * (w * t).cos();
            let zme = zm;

            let distrmnm = distance3dkm(xme, yme, zme, xtinit, ytinit, ztinit);
            let altmnm = ((xm * xm + ym * ym + zm * zm).sqrt() - a) / 3280.0;
            let axmguidg = (axmguid * axmguid + aymguid * aymguid + azmguid * azmguid).sqrt() / 32.2;

            array_t.push(t);
            array_distrtnm.push(distrtnm);
            array_altnm.push(altnm);
            array_distrmnm.push(distrmnm);
            array_altmnm.push(altmnm);
            array_axmguidg.push(axmguidg);
        }
    }

    // Final output point
    let xte = xt * (w * t).cos() - yt * (w * t).sin();
    let yte = xt * (w * t).sin() + yt * (w * t).cos();
    let zte = zt;
    let distrtnm = distance3dkm(xte, yte, zte, xtinit, ytinit, ztinit);
    let altnm = ((xt * xt + yt * yt + zt * zt).sqrt() - a) / 3280.0;
    let xme = xm * (w * t).cos() - ym * (w * t).sin();
    let yme = xm * (w * t).sin() + ym * (w * t).cos();
    let zme = zm;
    let distrmnm = distance3dkm(xme, yme, zme, xtinit, ytinit, ztinit);
    let altmnm = ((xm * xm + ym * ym + zm * zm).sqrt() - a) / 3280.0;
    let axmguidg = (axmguid * axmguid + aymguid * aymguid + azmguid * azmguid).sqrt() / 32.2;

    array_t.push(t);
    array_distrtnm.push(distrtnm);
    array_altnm.push(altnm);
    array_distrmnm.push(distrmnm);
    array_altmnm.push(altmnm);
    array_axmguidg.push(axmguidg);

    Results {
        time: array_t,
        dist_rt_nm: array_distrtnm,
        alt_nm: array_altnm,
        dist_rm_nm: array_distrmnm,
        alt_m_nm: array_altmnm,
        axm_guid_g: array_axmguidg,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c40l4_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.dist_rt_nm.clone(),
        results.alt_nm.clone(),
        results.dist_rm_nm.clone(),
        results.alt_m_nm.clone(),
        results.axm_guid_g.clone(),
    ])?;

    println!("C40L4: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c40l4_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
