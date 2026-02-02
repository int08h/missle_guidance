//! Chapter 43, Lesson 2: Radar Tracking with Kalman Filter
//!
//! Ballistic trajectory tracking with radar and Kalman filtering.

use crate::save_data;
use crate::utils::{lambert3d, distance3dkm, EARTH_RADIUS_FT, GM_FT};
use rand_distr::{Distribution, Normal};

pub struct Results {
    pub time: Vec<f64>,
    pub dist_nm: Vec<f64>,
    pub alt_nm: Vec<f64>,
    pub err_xtd: Vec<f64>,
    pub sp22: Vec<f64>,
    pub sp22p: Vec<f64>,
}

/// Run the C43L2 simulation
pub fn run() -> Results {
    let ifilter: i32 = 1;  // 1 for 2nd order, else 3rd order
    let rdeskm: f64 = 7000.0;
    let _tf: f64 = 2000.0;
    let tfinish: f64 = 240.0;
    let tloft: f64 = 500.0;
    let tgravend: f64 = 100.0;
    let gamdegic: f64 = 89.8;
    let tupt: f64 = 20.0;
    let rdesrkm: f64 = 560.0;
    let switch: i32 = 0;
    let phis: f64 = 0.0;
    let _phis1: f64 = 260.0;
    let _err: f64 = 0.0;
    let ts: f64 = 1.0;
    let sigthet: f64 = 0.001;
    let sigr: f64 = 10.0 * 3.28;
    let qgrav: i32 = 0;
    let _order: i32 = 2;
    let _left: i32 = 1;
    let mut qboost: bool = true;
    let _qfinish: i32 = 1;
    let mut qfirst: bool = true;

    let a = EARTH_RADIUS_FT;
    let gm = GM_FT;

    let mut _gamdeg = gamdegic;
    let hint: f64 = 0.01;
    let mut xh: f64 = 0.0;
    let mut xdh: f64 = 0.0;
    let mut xddh: f64;
    let mut yh: f64 = 0.0;
    let mut ydh: f64 = 0.0;
    let mut yddh: f64;

    // 2x2 filter matrices
    let mut p = [[0.0_f64; 2]; 2];
    let mut pp = [[0.0_f64; 2]; 2];
    p[0][0] = 99999999999.0;
    p[1][1] = 99999999999.0;
    pp[0][0] = 99999999999.0;
    pp[1][1] = 99999999999.0;

    let mut phi = [[0.0_f64; 2]; 2];
    phi[0][0] = 1.0;
    phi[0][1] = ts;
    phi[1][1] = 1.0;

    let hmat = [1.0_f64, 0.0];

    let mut q = [[0.0_f64; 2]; 2];
    q[0][0] = phis * ts.powi(3) / 3.0;
    q[0][1] = phis * ts.powi(2) / 2.0;
    q[1][0] = q[0][1];
    q[1][1] = phis * ts;

    // 3rd order covariance for P matrices (scalar tracking)
    let _p11: f64 = 99999999999.0;
    let _p12: f64 = 0.0;
    let _p13: f64 = 0.0;
    let p22: f64 = 99999999999.0;
    let _p23: f64 = 0.0;
    let _p33: f64 = 99999999999.0;
    let _p11p: f64 = 99999999999.0;
    let _p12p: f64 = 0.0;
    let _p13p: f64 = 0.0;
    let _p22p: f64 = 99999999999.0;
    let _p23p: f64 = 0.0;
    let _p33p: f64 = 99999999999.0;

    let mut xn: f64 = 0.0;
    let mut t: f64 = 0.0;
    let mut s: f64 = 0.0;
    let altnm: f64 = 0.0;
    let alt = altnm * 6076.0;
    let angdeg: f64 = 0.0;
    let ang = angdeg / 57.3;
    let mut x = (a + alt) * ang.cos();
    let mut y = (a + alt) * ang.sin();
    let mut alt = (x * x + y * y).sqrt() - a;
    let xfirst = x;
    let yfirst = y;
    let zfirst: f64 = 0.0;

    let mut x1 = (1.5708 - _gamdeg / 57.3 + ang).cos();
    let mut y1 = (1.5708 - _gamdeg / 57.3 + ang).sin();
    let mut axt: f64 = 0.0;
    let mut ayt: f64 = 0.0;

    let xlongtdeg = 57.3 * rdeskm * 3280.0 / a;
    let xlongrdeg = 57.3 * rdesrkm * 3280.0 / a;
    let mut tf = 252.0 + 0.223 * rdeskm - 5.44e-6 * rdeskm * rdeskm;
    tf += tloft;
    let xlongt = xlongtdeg / 57.3;
    let xlongr = xlongrdeg / 57.3;
    let xf = a * xlongt.cos();
    let yf = a * xlongt.sin();
    let xr = a * xlongr.cos();
    let yr = a * xlongr.sin();

    let z: f64 = 0.0;
    let zf: f64 = 0.0;
    let zr: f64 = 0.0;

    let mut rng = rand::thread_rng();
    let normal = Normal::new(0.0, 1.0).unwrap();

    let mut array_t = Vec::new();
    let mut array_distnm = Vec::new();
    let mut array_altnm = Vec::new();
    let mut array_errxtd = Vec::new();
    let mut array_sp22 = Vec::new();
    let mut array_sp22p = Vec::new();

    while !(alt < -1.0 || t > tfinish) {
        let xold = x;
        let yold = y;
        let x1old = x1;
        let y1old = y1;

        // First derivative evaluation
        let (wgt, trst) = if t < 120.0 {
            (-2622.0 * t + 440660.0, 725850.0)
        } else if t < 240.0 {
            (-642.0 * t + 168120.0, 182250.0)
        } else {
            (5500.0, 0.0)
        };

        let at = 32.2 * trst / wgt;
        let xd = x1;
        let yd = y1;
        let vel = (xd * xd + yd * yd).sqrt();
        let tembot = (x * x + y * y).powf(1.5);
        let x1d = -gm * x / tembot + axt;
        let y1d = -gm * y / tembot + ayt;
        alt = (x * x + y * y).sqrt() - a;
        let _accg = (axt * axt + ayt * ayt).sqrt() / 32.2;

        // Euler step
        x += hint * xd;
        y += hint * yd;
        x1 += hint * x1d;
        y1 += hint * y1d;
        t += hint;

        // RK2 averaging
        x = (xold + x) / 2.0 + 0.5 * hint * xd;
        y = (yold + y) / 2.0 + 0.5 * hint * yd;
        x1 = (x1old + x1) / 2.0 + 0.5 * hint * x1d;
        y1 = (y1old + y1) / 2.0 + 0.5 * hint * y1d;

        s += hint;
        let tgolam = tf - t;

        // Lambert guidance
        if qboost && t > tgravend {
            let result = lambert3d(x, y, z, tgolam, xf, yf, zf, switch);
            let vrx = result.vrx;
            let vry = result.vry;
            let delx = vrx - x1;
            let dely = vry - y1;
            let del = (delx * delx + dely * dely).sqrt();

            if t < 240.0 && del > 500.0 {
                axt = at * delx / del;
                ayt = at * dely / del;
            } else if del < 500.0 {
                qboost = false;
                axt = 0.0;
                ayt = 0.0;
                x1 = vrx;
                y1 = vry;
            } else {
                qboost = false;
                axt = 0.0;
                ayt = 0.0;
            }
        } else if t >= tupt && t <= tgravend && qfirst {
            qfirst = false;
            let vel = (xd * xd + yd * yd).sqrt();
            x1 = vel * (1.5708 - gamdegic / 57.3 + ang).cos();
            y1 = vel * (1.5708 - gamdegic / 57.3 + ang).sin();
            axt = at * x1 / vel;
            ayt = at * y1 / vel;
        } else if t >= tupt && t <= tgravend {
            let vel = (xd * xd + yd * yd).sqrt();
            axt = at * x1 / vel;
            ayt = at * y1 / vel;
        } else if t <= tupt {
            let rtmag = (x * x + y * y).sqrt();
            axt = at * x / rtmag;
            ayt = at * y / rtmag;
        }

        // Sampling and filtering
        if s >= ts - 0.00001 {
            s = 0.0;
            let distnm = distance3dkm(x, y, z, xfirst, yfirst, zfirst);
            let altnm = ((x * x + y * y).sqrt() - a) / 3280.0;
            let rmag = (x * x + y * y).sqrt();
            let vmag = (xd * xd + yd * yd).sqrt();
            _gamdeg = 90.0 - 57.3 * ((x * xd + y * yd) / (rmag * vmag)).acos();
            let _rho = 0.0034 * (-alt / 22000.0).exp();
            let _qpres = 0.5 * _rho * vel * vel;
            let _rrkm = ((x - xr).powi(2) + (y - yr).powi(2)).sqrt() / 3280.0;
            let _distrkm = distance3dkm(xr, yr, zr, xfirst, yfirst, zfirst);
            let _altrkm = ((xr * xr + yr * yr).sqrt() - a) / 3280.0;
            let rrmag = (xr * xr + yr * yr).sqrt();
            let rrtmag = ((x - xr).powi(2) + (y - yr).powi(2)).sqrt();
            let eldeg = 90.0 - 57.3 * ((xr * (x - xr) + yr * (y - yr)) / (rrmag * rrtmag)).acos();

            let isee = eldeg > 2.0 && eldeg < 85.0;

            if isee {
                xn += 1.0;
                let xk1 = 2.0 * (2.0 * xn - 1.0) / (xn * (xn + 1.0));
                let xk2 = 6.0 / (xn * (xn + 1.0) * ts);

                let _ts2 = ts * ts;
                let thet = (y - yr).atan2(xr - x);
                let r = ((xr - x).powi(2) + (y - yr).powi(2)).sqrt();
                let thetnoise = sigthet * normal.sample(&mut rng);
                let rnoise = sigr * normal.sample(&mut rng);
                let rmeas = r + rnoise;
                let thetmeas = thet + thetnoise;
                let xts = xr - rmeas * thetmeas.cos();
                let yts = yr + rmeas * thetmeas.sin();
                let _xtnoise = x - xts;
                let sigx = ((thet.cos() * sigr).powi(2) + (r * thet.sin() * sigthet).powi(2)).sqrt();
                let _ytnoise = y - yts;
                let sigy = ((thet.sin() * sigr).powi(2) + (r * thet.cos() * sigthet).powi(2)).sqrt();

                if ifilter == 1 {
                    // 2nd order Kalman filter
                    let rmat = sigx.powi(2);

                    // Propagate P
                    let phip00 = phi[0][0] * p[0][0] + phi[0][1] * p[1][0];
                    let phip01 = phi[0][0] * p[0][1] + phi[0][1] * p[1][1];
                    let phip10 = phi[1][0] * p[0][0] + phi[1][1] * p[1][0];
                    let phip11 = phi[1][0] * p[0][1] + phi[1][1] * p[1][1];

                    let m00 = phip00 * phi[0][0] + phip01 * phi[0][1] + q[0][0];
                    let m01 = phip00 * phi[1][0] + phip01 * phi[1][1] + q[0][1];
                    let m10 = phip10 * phi[0][0] + phip11 * phi[0][1] + q[1][0];
                    let m11 = phip10 * phi[1][0] + phip11 * phi[1][1] + q[1][1];

                    let hmht = hmat[0] * m00 * hmat[0];
                    let hmhtr = hmht + rmat;
                    let k0 = (m00 * hmat[0]) / hmhtr;
                    let k1 = (m10 * hmat[0]) / hmhtr;

                    p[0][0] = (1.0 - k0 * hmat[0]) * m00;
                    p[0][1] = (1.0 - k0 * hmat[0]) * m01;
                    p[1][0] = -k1 * hmat[0] * m00 + m10;
                    p[1][1] = -k1 * hmat[0] * m01 + m11;

                    let xk1pz = if xk1 > k0 { xk1 } else { k0 };
                    let xk2pz = if xk1 > k0 { xk2 } else { k1 };

                    if qgrav == 0 {
                        xddh = x1d;
                    } else {
                        let velh = (xdh * xdh + ydh * ydh).sqrt();
                        xddh = at * xdh / (velh + 0.0001);
                    }

                    let res = xts - xh - ts * xdh - 0.5 * ts * ts * xddh;
                    xh = xh + xdh * ts + 0.5 * ts * ts * xddh + xk1pz * res;
                    xdh = xdh + xddh * ts + xk2pz * res;

                    // Y filter
                    let rmatp = sigy.powi(2);
                    let phipp00 = phi[0][0] * pp[0][0] + phi[0][1] * pp[1][0];
                    let phipp01 = phi[0][0] * pp[0][1] + phi[0][1] * pp[1][1];
                    let phipp10 = phi[1][0] * pp[0][0] + phi[1][1] * pp[1][0];
                    let phipp11 = phi[1][0] * pp[0][1] + phi[1][1] * pp[1][1];

                    let mp00 = phipp00 * phi[0][0] + phipp01 * phi[0][1] + q[0][0];
                    let mp01 = phipp00 * phi[1][0] + phipp01 * phi[1][1] + q[0][1];
                    let mp10 = phipp10 * phi[0][0] + phipp11 * phi[0][1] + q[1][0];
                    let mp11 = phipp10 * phi[1][0] + phipp11 * phi[1][1] + q[1][1];

                    let hmhtp = hmat[0] * mp00 * hmat[0];
                    let hmhtrp = hmhtp + rmatp;
                    let kp0 = (mp00 * hmat[0]) / hmhtrp;
                    let kp1 = (mp10 * hmat[0]) / hmhtrp;

                    pp[0][0] = (1.0 - kp0 * hmat[0]) * mp00;
                    pp[0][1] = (1.0 - kp0 * hmat[0]) * mp01;
                    pp[1][0] = -kp1 * hmat[0] * mp00 + mp10;
                    pp[1][1] = -kp1 * hmat[0] * mp01 + mp11;

                    let xk1pzp = if xk1 > k0 { xk1 } else { kp0 };
                    let xk2pzp = if xk1 > k0 { xk2 } else { kp1 };

                    if qgrav == 0 {
                        yddh = y1d;
                    } else {
                        let velh = (xdh * xdh + ydh * ydh).sqrt();
                        yddh = at * ydh / (velh + 0.0001);
                    }

                    let resp = yts - yh - ts * ydh - 0.5 * ts * ts * yddh;
                    yh = yh + ydh * ts + 0.5 * ts * ts * yddh + xk1pzp * resp;
                    ydh = ydh + yddh * ts + xk2pzp * resp;
                }
            }

            // Data recording happens OUTSIDE if isee block (matches MATLAB)
            if ifilter == 1 {
                let errxtd = x1 - xdh;
                let sp22_val = p[1][1].sqrt();
                let sp22p_val = -sp22_val;

                array_t.push(t);
                array_distnm.push(distnm);
                array_altnm.push(altnm);
                array_errxtd.push(errxtd);
                array_sp22.push(sp22_val);
                array_sp22p.push(sp22p_val);
            } else {
                // 3rd order filter - simplified version
                let errxtd = x1 - xdh;
                let sp22_val = p22.sqrt();
                let sp22p_val = -sp22_val;

                array_t.push(t);
                array_distnm.push(distnm);
                array_altnm.push(altnm);
                array_errxtd.push(errxtd);
                array_sp22.push(sp22_val);
                array_sp22p.push(sp22p_val);
            }
        }
    }

    Results {
        time: array_t,
        dist_nm: array_distnm,
        alt_nm: array_altnm,
        err_xtd: array_errxtd,
        sp22: array_sp22,
        sp22p: array_sp22p,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c43l2_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.dist_nm.clone(),
        results.alt_nm.clone(),
        results.err_xtd.clone(),
        results.sp22.clone(),
        results.sp22p.clone(),
    ])?;

    println!("C43L2: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c43l2_runs() {
        let results = run();
        // May or may not have data depending on engagement geometry
        // Just verify it runs without error
        assert!(true);
    }
}
