//! Chapter 45, Lesson 3: TMD Intercept Simulation
//!
//! Theater missile defense: intercept simulation with guidance.

use crate::save_data;
use crate::utils::{lambert3d, distance3dkm, predict45, EARTH_RADIUS_FT, GM_FT};

pub struct Results {
    pub time: Vec<f64>,
    pub dist_rt_km: Vec<f64>,
    pub alt_t_km: Vec<f64>,
    pub dist_rm_km: Vec<f64>,
    pub alt_m_km: Vec<f64>,
    pub atplos_g: Vec<f64>,
    pub xncplos_g: Vec<f64>,
}

/// Run the C45L3 simulation
pub fn run() -> Results {
    let tlaunch: f64 = 90.0;
    let mut tf: f64 = 170.0;
    let ts: f64 = 1.0;
    let xlongmdegickm: f64 = 400.0;
    let rdeskm: f64 = 2000.0;
    let gamdeg: f64 = 89.99;
    let tupt: f64 = 15.0;
    let tguid: f64 = 110.0;
    let xnclim: f64 = 322.0;
    let xnp: f64 = 3.0;
    let qperfect: bool = true;
    let tloft: f64 = 200.0;
    let altmkmic: f64 = 15.0;
    let qtaylor: bool = true;
    let qguid: bool = true;
    let itgt: i32 = 1;
    let switch1: i32 = 0;
    let switchm: i32 = 0;
    let deltf: f64 = 0.0;
    let _qfix: bool = true;

    let tpz = if itgt == 1 { 180.0 } else { 240.0 };

    let a = EARTH_RADIUS_FT;
    let gm = GM_FT;

    let altm = altmkmic * 3280.0;
    let mut tftot = 252.0 + 0.223 * rdeskm - 5.44e-6 * rdeskm * rdeskm;
    tftot += tloft;

    let xlongfdeg = 57.3 * rdeskm * 3280.0 / a;
    let xlongmdegic = xlongmdegickm / 111.0;
    let xlongmdeg = xlongmdegic;
    let xlongtdeg: f64 = 0.0;

    let xlongf = xlongfdeg / 57.3;
    let xlongt = xlongtdeg / 57.3;
    let xlongm = xlongmdeg / 57.3;

    let xf = a * xlongf.cos();
    let yf = a * xlongf.sin();
    let zf: f64 = 0.0;

    let mut xt = a * xlongt.cos();
    let mut yt = a * xlongt.sin();
    let zt: f64 = 0.0;

    let mut xm = (a + altm) * xlongm.cos();
    let mut ym = (a + altm) * xlongm.sin();
    let zm: f64 = 0.0;
    let xfirst = xt;
    let yfirst = yt;
    let zfirst = zt;

    let mut xtd = (1.5708 - gamdeg / 57.3).cos();
    let mut ytd = (1.5708 - gamdeg / 57.3).sin();
    let mut xmd: f64 = 0.0;
    let mut ymd: f64 = 0.0;

    let mut rtm1 = xt - xm;
    let mut rtm2 = yt - ym;
    let mut rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
    let vtm1 = xtd - xmd;
    let vtm2 = ytd - ymd;
    let mut vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;

    let mut delv: f64 = 0.0;
    let mut axmguid: f64;
    let mut aymguid: f64;
    let mut zem1: f64;
    let mut zem2: f64;
    let mut axt: f64 = 0.0;
    let mut ayt: f64 = 0.0;

    // Find exact target location at intercept time
    let result_pred = predict45(0.0, xt, yt, xtd, ytd, tf, tftot, tupt, xf, yf, itgt);
    let xtfact = result_pred.xtf;
    let ytfact = result_pred.ytf;

    // Find required missile velocity
    let tgolam = tf - tlaunch;
    let result_lam = lambert3d(xm, ym, zm, tgolam, xtfact, ytfact, 0.0, switchm);
    let _vmxrqd = result_lam.vrx;
    let _vmyrqd = result_lam.vry;

    tf += deltf;

    let mut t: f64 = 0.0;
    let mut s: f64 = 0.0;
    let mut h: f64;
    let mut qboost: bool = true;
    let mut qfirst: bool = true;

    let mut array_t = Vec::new();
    let mut array_distrtkm = Vec::new();
    let mut array_alttkm = Vec::new();
    let mut array_distrmkm = Vec::new();
    let mut array_altmkm = Vec::new();
    let mut array_atplosg = Vec::new();
    let mut array_xncplosg = Vec::new();

    while !((t > tf - 10.0) && vc < 0.0) {
        if rtm < 1000.0 {
            h = 0.00001;
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
        let delvold = delv;

        // First derivative evaluation
        let (_wgt, _trst) = if itgt == 1 {
            if t < 180.0 {
                (-212.0 * t + 44000.0, 54100.0)
            } else {
                (3300.0, 0.0)
            }
        } else if t < 120.0 {
            (-2622.0 * t + 440660.0, 725850.0)
        } else if t < 240.0 {
            (-642.0 * t + 168120.0, 182250.0)
        } else {
            (5500.0, 0.0)
        };

        let tempbott = (xt * xt + yt * yt).powf(1.5);
        let mut xtdd = -gm * xt / tempbott + axt;
        let mut ytdd = -gm * yt / tempbott + ayt;

        rtm1 = xt - xm;
        rtm2 = yt - ym;
        let vtm1 = xtd - xmd;
        let vtm2 = ytd - ymd;
        rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
        vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;
        let tgo = rtm / vc;

        // ZEM guidance
        if t > tguid {
            let tempbotm = (xm * xm + ym * ym).powf(1.5);
            let xmddgrav = -gm * xm / tempbotm;
            let ymddgrav = -gm * ym / tempbotm;
            zem1 = rtm1 + vtm1 * tgo + 0.5 * (xtdd - xmddgrav) * tgo * tgo;
            zem2 = rtm2 + vtm2 * tgo + 0.5 * (ytdd - ymddgrav) * tgo * tgo;
            let zemdotrtm = (zem1 * rtm1 + zem2 * rtm2) / rtm;
            let zemper1 = zem1 - zemdotrtm * rtm1 / rtm;
            let zemper2 = zem2 - zemdotrtm * rtm2 / rtm;

            axmguid = xnp * zemper1 / (tgo * tgo);
            aymguid = xnp * zemper2 / (tgo * tgo);

            let lim = if qguid { xnclim } else { 0.0 };
            axmguid = axmguid.clamp(-lim, lim);
            aymguid = aymguid.clamp(-lim, lim);
        } else {
            axmguid = 0.0;
            aymguid = 0.0;
        }

        let (xmdd, ymdd) = if t > tlaunch {
            let tempbotm = (xm * xm + ym * ym).powf(1.5);
            (-gm * xm / tempbotm + axmguid, -gm * ym / tempbotm + aymguid)
        } else {
            (0.0, 0.0)
        };

        let accnew = (axmguid * axmguid + aymguid * aymguid).sqrt();
        let delvd = accnew;
        let _altmkm = ((xm * xm + ym * ym).sqrt() - a) / 3280.0;

        // Euler step
        xt += h * xtd;
        yt += h * ytd;
        xtd += h * xtdd;
        ytd += h * ytdd;
        xm += h * xmd;
        ym += h * ymd;
        xmd += h * xmdd;
        ymd += h * ymdd;
        delv += h * delvd;
        t += h;

        // Second derivative evaluation for RK2
        let (wgt2, trst2) = if itgt == 1 {
            if t < 180.0 {
                (-212.0 * t + 44000.0, 54100.0)
            } else {
                (3300.0, 0.0)
            }
        } else if t < 120.0 {
            (-2622.0 * t + 440660.0, 725850.0)
        } else if t < 240.0 {
            (-642.0 * t + 168120.0, 182250.0)
        } else {
            (5500.0, 0.0)
        };

        let atp = 32.2 * trst2 / wgt2;
        let tempbott2 = (xt * xt + yt * yt).powf(1.5);
        xtdd = -gm * xt / tempbott2 + axt;
        ytdd = -gm * yt / tempbott2 + ayt;

        rtm1 = xt - xm;
        rtm2 = yt - ym;
        let vtm1 = xtd - xmd;
        let vtm2 = ytd - ymd;
        rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
        vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;
        let tgo2 = rtm / vc;

        if t > tguid {
            let tempbotm = (xm * xm + ym * ym).powf(1.5);
            let xmddgrav = -gm * xm / tempbotm;
            let ymddgrav = -gm * ym / tempbotm;
            zem1 = rtm1 + vtm1 * tgo2 + 0.5 * (xtdd - xmddgrav) * tgo2 * tgo2;
            zem2 = rtm2 + vtm2 * tgo2 + 0.5 * (ytdd - ymddgrav) * tgo2 * tgo2;
            let zemdotrtm = (zem1 * rtm1 + zem2 * rtm2) / rtm;
            let zemper1 = zem1 - zemdotrtm * rtm1 / rtm;
            let zemper2 = zem2 - zemdotrtm * rtm2 / rtm;

            axmguid = xnp * zemper1 / (tgo2 * tgo2);
            aymguid = xnp * zemper2 / (tgo2 * tgo2);

            let lim = if qguid { xnclim } else { 0.0 };
            axmguid = axmguid.clamp(-lim, lim);
            aymguid = aymguid.clamp(-lim, lim);
        } else {
            axmguid = 0.0;
            aymguid = 0.0;
        }

        let (xmdd2, ymdd2) = if t > tlaunch {
            let tempbotm = (xm * xm + ym * ym).powf(1.5);
            (-gm * xm / tempbotm + axmguid, -gm * ym / tempbotm + aymguid)
        } else {
            (0.0, 0.0)
        };

        let accnew2 = (axmguid * axmguid + aymguid * aymguid).sqrt();
        let delvd2 = accnew2;

        // RK2 averaging
        xt = 0.5 * (xtold + xt + h * xtd);
        yt = 0.5 * (ytold + yt + h * ytd);
        xtd = 0.5 * (xtdold + xtd + h * xtdd);
        ytd = 0.5 * (ytdold + ytd + h * ytdd);
        xm = 0.5 * (xmold + xm + h * xmd);
        ym = 0.5 * (ymold + ym + h * ymd);
        xmd = 0.5 * (xmdold + xmd + h * xmdd2);
        ymd = 0.5 * (ymdold + ymd + h * ymdd2);
        delv = 0.5 * (delvold + delv + h * delvd2);

        s += h;

        // Lambert guidance for target
        if qboost {
            let tgolam = tftot - t;
            let result = lambert3d(xt, yt, zt, tgolam, xf, yf, zf, switch1);
            let vtx = result.vrx;
            let vty = result.vry;

            let delvxt = vtx - xtd;
            let delvyt = vty - ytd;
            let _velt = (xtd * xtd + ytd * ytd).sqrt();
            let delvelt = (delvxt * delvxt + delvyt * delvyt).sqrt();

            if t < tpz && delvelt > 500.0 {
                axt = atp * delvxt / delvelt;
                ayt = atp * delvyt / delvelt;
            } else if delvelt < 500.0 {
                qboost = false;
                axt = 0.0;
                ayt = 0.0;
                xtd = vtx;
                ytd = vty;
            } else {
                qboost = false;
                axt = 0.0;
                ayt = 0.0;
            }
        }

        if t < tupt {
            let rtmag = (xt * xt + yt * yt).sqrt();
            axt = atp * xt / rtmag;
            ayt = atp * yt / rtmag;
        }

        // Predicted intercept point
        let (xtf, ytf) = if t >= tlaunch {
            if qperfect {
                (xtfact, ytfact)
            } else if qtaylor {
                let tgom = tf - t;
                (xt + xtd * tgom + 0.5 * xtdd * tgom * tgom,
                 yt + ytd * tgom + 0.5 * ytdd * tgom * tgom)
            } else {
                (xtfact, ytfact)
            }
        } else {
            (xtfact, ytfact)
        };

        // Launch missile
        if t >= tlaunch && qfirst {
            qfirst = false;
            let tgopz = tf - tlaunch;
            let result = lambert3d(xm, ym, zm, tgopz, xtf, ytf, 0.0, switchm);
            xmd = result.vrx;
            ymd = result.vry;
        }

        if s >= ts - 0.0001 {
            s = 0.0;
            let alttkm = ((xt * xt + yt * yt).sqrt() - a) / 3280.0;
            let distrtkm = distance3dkm(xt, yt, zt, xfirst, yfirst, zfirst);
            let altmkm_val = ((xm * xm + ym * ym).sqrt() - a) / 3280.0;
            let distrmkm = distance3dkm(xm, ym, zm, xfirst, yfirst, zfirst);

            let xlam = rtm2.atan2(rtm1);
            let atplos = -xtdd * xlam.sin() + ytdd * xlam.cos();
            let atplosg = atplos / 32.2;
            let xncplosg = (-axmguid * xlam.sin() + aymguid * xlam.cos()) / 32.2;

            array_t.push(t);
            array_distrtkm.push(distrtkm);
            array_alttkm.push(alttkm);
            array_distrmkm.push(distrmkm);
            array_altmkm.push(altmkm_val);
            array_atplosg.push(atplosg);
            array_xncplosg.push(xncplosg);
        }
    }

    Results {
        time: array_t,
        dist_rt_km: array_distrtkm,
        alt_t_km: array_alttkm,
        dist_rm_km: array_distrmkm,
        alt_m_km: array_altmkm,
        atplos_g: array_atplosg,
        xncplos_g: array_xncplosg,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c45l3_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.dist_rt_km.clone(),
        results.alt_t_km.clone(),
        results.dist_rm_km.clone(),
        results.alt_m_km.clone(),
        results.atplos_g.clone(),
        results.xncplos_g.clone(),
    ])?;

    println!("C45L3: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c45l3_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
