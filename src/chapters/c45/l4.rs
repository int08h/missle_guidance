//! Chapter 45, Lesson 4: TMD Monte Carlo Simulation
//!
//! Theater missile defense: Monte Carlo simulation with tracking noise.

use crate::save_data;
use crate::utils::{lambert3d, predict45, EARTH_RADIUS_FT, GM_FT};
use rand_distr::{Distribution, Normal};

pub struct Results {
    pub run: Vec<f64>,
    pub rtm: Vec<f64>,
    pub delv: Vec<f64>,
}

/// Run the C45L4 simulation
pub fn run() -> Results {
    let phis: f64 = 576.0;
    let tlaunch: f64 = 80.0;
    let tf_init: f64 = 230.0;
    let ts: f64 = 0.1;
    let xlongmdegickm: f64 = 400.0;
    let xlongs2degkm: f64 = 500.0;
    let rdeskm: f64 = 10000.0;
    let gamdeg: f64 = 89.99;
    let tupt: f64 = 20.0;
    let tguid: f64 = 100.0;
    let xnclim: f64 = 322.0;
    let xnp: f64 = 3.0;
    let qperfect: bool = false;
    let tloft: f64 = 500.0;
    let altmkmic: f64 = 15.0;
    let qtaylor: bool = true;
    let qguid: bool = true;
    let itgt: i32 = 2;  // ICBM
    let deltf: f64 = 0.0;
    let qfix: bool = false;
    let sigthet1: f64 = 0.00005;
    let thom: f64 = 10.0;
    let _bias1: f64 = 0.0;
    let run_count: usize = 50;
    let vmrqdkmic: f64 = 4.0;

    let tpz = if itgt == 1 { 180.0 } else { 240.0 };
    let sigthet2 = sigthet1;
    let xlongs1degkm = xlongmdegickm;

    let a = EARTH_RADIUS_FT;
    let gm = GM_FT;

    let altm = altmkmic * 3280.0;
    let mut tftot = 252.0 + 0.223 * rdeskm - 5.44e-6 * rdeskm * rdeskm;
    tftot += tloft;

    let mut rng = rand::thread_rng();
    let normal = Normal::new(0.0, 1.0).unwrap();

    let mut array_jj = Vec::new();
    let mut array_rtm = Vec::new();
    let mut array_delv = Vec::new();

    for jj in 0..run_count {
        let switchm: i32 = 0;
        let switch1: i32 = 0;

        let xlongfdeg = 57.3 * rdeskm * 3280.0 / a;
        let xlongmdegic = xlongmdegickm / 111.0;
        let xlongmdeg = xlongmdegic;
        let xlongs1deg = xlongs1degkm / 111.0;
        let xlongs2deg = xlongs2degkm / 111.0;
        let xlongtdeg: f64 = 0.0;

        let xlongf = xlongfdeg / 57.3;
        let xlongt = xlongtdeg / 57.3;
        let xlongm = xlongmdeg / 57.3;
        let xlongs1 = xlongs1deg / 57.3;
        let xlongs2 = xlongs2deg / 57.3;

        let xf = a * xlongf.cos();
        let yf = a * xlongf.sin();
        let zf: f64 = 0.0;

        let mut xt = a * xlongt.cos();
        let mut yt = a * xlongt.sin();
        let zt: f64 = 0.0;

        let mut xm = (a + altm) * xlongm.cos();
        let mut ym = (a + altm) * xlongm.sin();
        let zm: f64 = 0.0;
        let xs1 = (a + altm) * xlongs1.cos();
        let ys1 = (a + altm) * xlongs1.sin();
        let xs2 = (a + altm) * xlongs2.cos();
        let ys2 = (a + altm) * xlongs2.sin();

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
        let mut axmguid: f64 = 0.0;
        let mut aymguid: f64 = 0.0;
        let mut axt: f64 = 0.0;
        let mut ayt: f64 = 0.0;

        // Find intercept time
        let mut tf = tf_init;
        if !qfix {
            let mut tf_search = tlaunch + 30.0;
            while tf_search <= tpz - 10.0 {
                let result_pred = predict45(0.0, xt, yt, xtd, ytd, tf_search, tftot, tupt, xf, yf, itgt);
                let xtfact = result_pred.xtf;
                let ytfact = result_pred.ytf;
                let tgolam = tf_search - tlaunch;
                let result = lambert3d(xm, ym, zm, tgolam, xtfact, ytfact, 0.0, switchm);
                let vmxrqd = result.vrx;
                let vmyrqd = result.vry;
                let vmrqdkm = (vmxrqd * vmxrqd + vmyrqd * vmyrqd).sqrt() / 3280.0;
                if vmrqdkm < vmrqdkmic {
                    tf = tf_search;
                    break;
                }
                tf_search += 10.0;
            }
        }

        tf += deltf;

        // Initialize Kalman filter state
        let mut xh: f64 = 0.0;
        let mut xdh: f64 = 0.0;
        let mut xddh: f64 = 0.0;
        let mut yh: f64 = 0.0;
        let mut ydh: f64 = 0.0;
        let mut yddh: f64 = 0.0;

        // Kalman filter covariance matrices (3x3)
        let mut p = [[0.0_f64; 3]; 3];
        let mut pp = [[0.0_f64; 3]; 3];
        p[0][0] = 99999999999.0;
        p[1][1] = 99999999999.0;
        p[2][2] = 99999999999.0;
        pp[0][0] = 99999999999.0;
        pp[1][1] = 99999999999.0;
        pp[2][2] = 99999999999.0;

        let mut phi = [[0.0_f64; 3]; 3];
        phi[0][0] = 1.0;
        phi[0][1] = ts;
        phi[0][2] = 0.5 * ts * ts;
        phi[1][1] = 1.0;
        phi[1][2] = ts;
        phi[2][2] = 1.0;

        let hmat = [1.0_f64, 0.0, 0.0];

        let ts2 = ts * ts;
        let ts3 = ts2 * ts;
        let ts4 = ts3 * ts;
        let ts5 = ts4 * ts;

        let mut q = [[0.0_f64; 3]; 3];
        q[0][0] = phis * ts5 / 20.0;
        q[0][1] = phis * ts4 / 8.0;
        q[0][2] = phis * ts3 / 6.0;
        q[1][0] = q[0][1];
        q[1][1] = phis * ts3 / 3.0;
        q[1][2] = phis * ts2 / 2.0;
        q[2][0] = q[0][2];
        q[2][1] = q[1][2];
        q[2][2] = phis * ts;

        let mut xn: f64 = 0.0;

        // Find predicted intercept point
        let result_pred = predict45(0.0, xt, yt, xtd, ytd, tf, tftot, tupt, xf, yf, itgt);
        let xtfact = result_pred.xtf;
        let ytfact = result_pred.ytf;

        let mut t: f64 = 0.0;
        let mut s: f64 = 0.0;
        let mut h: f64;
        let mut qboost: bool = true;
        let mut qfirst: bool = true;
        let mut xtdd: f64;
        let mut ytdd: f64;

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
            let (wgt, trst) = if itgt == 1 {
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

            let _atp = 32.2 * trst / wgt;
            let tempbott = (xt * xt + yt * yt).powf(1.5);
            xtdd = -gm * xt / tempbott + axt;
            ytdd = -gm * yt / tempbott + ayt;

            rtm1 = xt - xm;
            rtm2 = yt - ym;
            let vtm1 = xtd - xmd;
            let vtm2 = ytd - ymd;
            rtm = (rtm1 * rtm1 + rtm2 * rtm2).sqrt();
            vc = -(rtm1 * vtm1 + rtm2 * vtm2) / rtm;
            let tgo = rtm / vc;

            // ZEM guidance with homing condition
            if t > tguid && tgo < thom {
                let tempbotm = (xm * xm + ym * ym).powf(1.5);
                let xmddgrav = -gm * xm / tempbotm;
                let ymddgrav = -gm * ym / tempbotm;
                let zem1 = rtm1 + vtm1 * tgo + 0.5 * (xtdd - xmddgrav) * tgo * tgo;
                let zem2 = rtm2 + vtm2 * tgo + 0.5 * (ytdd - ymddgrav) * tgo * tgo;
                let zemdotrtm = (zem1 * rtm1 + zem2 * rtm2) / rtm;
                let zemper1 = zem1 - zemdotrtm * rtm1 / rtm;
                let zemper2 = zem2 - zemdotrtm * rtm2 / rtm;

                axmguid = xnp * zemper1 / (tgo * tgo);
                aymguid = xnp * zemper2 / (tgo * tgo);

                let lim = if qguid { xnclim } else { 0.0 };
                axmguid = axmguid.clamp(-lim, lim);
                aymguid = aymguid.clamp(-lim, lim);
            }

            if t <= tguid {
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

            if t > tguid && tgo2 < thom {
                let tempbotm = (xm * xm + ym * ym).powf(1.5);
                let xmddgrav = -gm * xm / tempbotm;
                let ymddgrav = -gm * ym / tempbotm;
                let zem1 = rtm1 + vtm1 * tgo2 + 0.5 * (xtdd - xmddgrav) * tgo2 * tgo2;
                let zem2 = rtm2 + vtm2 * tgo2 + 0.5 * (ytdd - ymddgrav) * tgo2 * tgo2;
                let zemdotrtm = (zem1 * rtm1 + zem2 * rtm2) / rtm;
                let zemper1 = zem1 - zemdotrtm * rtm1 / rtm;
                let zemper2 = zem2 - zemdotrtm * rtm2 / rtm;

                axmguid = xnp * zemper1 / (tgo2 * tgo2);
                aymguid = xnp * zemper2 / (tgo2 * tgo2);

                let lim = if qguid { xnclim } else { 0.0 };
                axmguid = axmguid.clamp(-lim, lim);
                aymguid = aymguid.clamp(-lim, lim);
            }

            if t <= tguid {
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

            // Kalman filter update at sample times
            if s >= ts - 0.0001 {
                s = 0.0;

                // Triangulation measurements
                let thet1 = (ys1 - yt).atan2(xs1 - xt);
                let thet2 = (ys2 - yt).atan2(xs2 - xt);
                let thet1noise = sigthet1 * normal.sample(&mut rng);
                let thet2noise = sigthet2 * normal.sample(&mut rng);
                let thet1s = thet1 + thet1noise;
                let thet2s = thet2 + thet2noise;

                // Triangulation solution
                let top1 = xs2 * thet2s.tan() - xs1 * thet1s.tan() + ys1 - ys2;
                let xts = top1 / (thet2s.tan() - thet1s.tan());
                let top2 = xs2 * thet2s.tan() * thet1s.tan() - xs1 * thet1s.tan() * thet1s.tan()
                    + thet1s.tan() * (ys1 - ys2);
                let yts = ys1 - xs1 * thet1s.tan() + top2 / (thet2s.tan() - thet1s.tan());

                // Compute measurement noise
                let dxdt1 = (thet2.tan() * (xs2 - xs1) + ys1 - ys2)
                    / (thet1.cos() * (thet2.tan() - thet1.tan())).powi(2);
                let dxdt2 = (thet1.tan() * (xs1 - xs2) + ys2 - ys1)
                    / (thet2.cos() * (thet2.tan() - thet1.tan())).powi(2);
                let sigx = ((dxdt1 * sigthet1).powi(2) + (dxdt2 * sigthet2).powi(2)).sqrt();

                let dydt1 = -xs1 / (thet1.cos() * thet1.cos())
                    + (xs2 * thet2.tan() * thet2.tan() - 2.0 * xs1 * thet1.tan() * thet2.tan()
                        + (ys1 - ys2) * thet2.tan() + xs1 * thet1.tan() * thet1.tan())
                        / (thet1.cos() * (thet2.tan() - thet1.tan())).powi(2);
                let dydt2 = (thet1.tan() * thet1.tan() * (xs1 - xs2) - (ys1 - ys2) * thet1.tan())
                    / (thet2.cos() * (thet2.tan() - thet1.tan())).powi(2);
                let sigy = ((dydt1 * sigthet1).powi(2) + (dydt2 * sigthet2).powi(2)).sqrt();

                xn += 1.0;
                let xk1 = 3.0 * (3.0 * xn * xn - 3.0 * xn + 2.0) / (xn * (xn + 1.0) * (xn + 2.0));
                let xk2 = 18.0 * (2.0 * xn - 1.0) / (xn * (xn + 1.0) * (xn + 2.0) * ts);
                let xk3 = 60.0 / (xn * (xn + 1.0) * (xn + 2.0) * ts * ts);

                // Kalman filter for X
                let rmat = sigx * sigx;

                // Propagate P: M = PHI * P * PHI' + Q
                let mut m = [[0.0_f64; 3]; 3];
                for i in 0..3 {
                    for j in 0..3 {
                        let mut sum = 0.0;
                        for k in 0..3 {
                            let mut phip = 0.0;
                            for l in 0..3 {
                                phip += phi[i][l] * p[l][k];
                            }
                            sum += phip * phi[j][k];
                        }
                        m[i][j] = sum + q[i][j];
                    }
                }

                // K = M * H' / (H * M * H' + R)
                let hmht = m[0][0] * hmat[0] * hmat[0];
                let hmhtr = hmht + rmat;
                let k = [m[0][0] * hmat[0] / hmhtr, m[1][0] * hmat[0] / hmhtr, m[2][0] * hmat[0] / hmhtr];

                // P = (I - K*H) * M
                for i in 0..3 {
                    for j in 0..3 {
                        p[i][j] = m[i][j] - k[i] * hmat[0] * m[0][j];
                    }
                }

                let (xk1pz, xk2pz, xk3pz) = if xn < 10.0 {
                    (xk1, xk2, xk3)
                } else {
                    (k[0], k[1], k[2])
                };

                let res = xts - xh - ts * xdh - 0.5 * ts * ts * xddh;
                xh = xh + xdh * ts + 0.5 * ts * ts * xddh + xk1pz * res;
                xdh = xdh + xddh * ts + xk2pz * res;
                xddh = xddh + xk3pz * res;

                // Kalman filter for Y
                let rmatp = sigy * sigy;

                let mut mp = [[0.0_f64; 3]; 3];
                for i in 0..3 {
                    for j in 0..3 {
                        let mut sum = 0.0;
                        for k in 0..3 {
                            let mut phip = 0.0;
                            for l in 0..3 {
                                phip += phi[i][l] * pp[l][k];
                            }
                            sum += phip * phi[j][k];
                        }
                        mp[i][j] = sum + q[i][j];
                    }
                }

                let hmhtp = mp[0][0] * hmat[0] * hmat[0];
                let hmhtrp = hmhtp + rmatp;
                let kp = [mp[0][0] * hmat[0] / hmhtrp, mp[1][0] * hmat[0] / hmhtrp, mp[2][0] * hmat[0] / hmhtrp];

                for i in 0..3 {
                    for j in 0..3 {
                        pp[i][j] = mp[i][j] - kp[i] * hmat[0] * mp[0][j];
                    }
                }

                let (xk1pzp, xk2pzp, xk3pzp) = if xn < 10.0 {
                    (xk1, xk2, xk3)
                } else {
                    (kp[0], kp[1], kp[2])
                };

                let resp = yts - yh - ts * ydh - 0.5 * ts * ts * yddh;
                yh = yh + ydh * ts + 0.5 * ts * ts * yddh + xk1pzp * resp;
                ydh = ydh + yddh * ts + xk2pzp * resp;
                yddh = yddh + xk3pzp * resp;

                // Use filtered state for guidance when TGO > THOM
                if t > tguid && tgo > thom {
                    let rtm1h = xh - xm;
                    let rtm2h = yh - ym;
                    let rtmh = (rtm1h * rtm1h + rtm2h * rtm2h).sqrt();
                    let vtm1h = xdh - xmd;
                    let vtm2h = ydh - ymd;
                    let vch = -(rtm1h * vtm1h + rtm2h * vtm2h) / rtmh;
                    let tgoh = rtmh / vch;
                    let tempbotm = (xm * xm + ym * ym).powf(1.5);
                    let xmddgrav = -gm * xm / tempbotm;
                    let ymddgrav = -gm * ym / tempbotm;
                    let zem1h = rtm1h + vtm1h * tgoh + 0.5 * (xddh - xmddgrav) * tgoh * tgoh;
                    let zem2h = rtm2h + vtm2h * tgoh + 0.5 * (yddh - ymddgrav) * tgoh * tgoh;
                    let zemdotrtmh = (zem1h * rtm1h + zem2h * rtm2h) / rtmh;
                    let zemper1h = zem1h - zemdotrtmh * rtm1h / rtmh;
                    let zemper2h = zem2h - zemdotrtmh * rtm2h / rtmh;

                    axmguid = xnp * zemper1h / (tgoh * tgoh);
                    aymguid = xnp * zemper2h / (tgoh * tgoh);

                    let lim = if qguid { xnclim } else { 0.0 };
                    axmguid = axmguid.clamp(-lim, lim);
                    aymguid = aymguid.clamp(-lim, lim);
                }
            }
        }

        array_jj.push((jj + 1) as f64);
        array_rtm.push(rtm);
        array_delv.push(delv / 3.28);
    }

    Results {
        run: array_jj,
        rtm: array_rtm,
        delv: array_delv,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c45l4_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.run.clone(),
        results.rtm.clone(),
        results.delv.clone(),
    ])?;

    println!("C45L4: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c45l4_runs() {
        let results = run();
        assert!(!results.run.is_empty());
    }
}
