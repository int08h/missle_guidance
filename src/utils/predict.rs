//! Trajectory prediction functions
//!
//! Implements predictpz, predict45, and predictg for ballistic and guided trajectory prediction

use super::constants::{EARTH_RADIUS_FT, GM_FT, G_ACCEL};
use super::lambert3d::lambert3d;

/// Result from pure ballistic prediction
#[derive(Debug, Clone, Copy)]
pub struct PredictResult {
    pub xtf: f64,
    pub ytf: f64,
}

/// Pure gravity (ballistic) trajectory prediction (predictpz.m)
///
/// Propagates position and velocity forward in time under gravity only.
///
/// # Arguments
/// * `tf` - Time to propagate
/// * `x`, `y` - Initial position
/// * `x1`, `y1` - Initial velocity
///
/// # Returns
/// Final position (xtf, ytf)
pub fn predict_pz(tf: f64, x: f64, y: f64, x1: f64, y1: f64) -> PredictResult {
    let h = 0.01;
    let gm = GM_FT;

    let mut t = 0.0;
    let mut x = x;
    let mut y = y;
    let mut x1 = x1;
    let mut y1 = y1;

    while t <= tf - 0.00001 {
        let xold = x;
        let yold = y;
        let x1old = x1;
        let y1old = y1;

        // First Euler step
        let tembot = (x * x + y * y).powf(1.5);
        let x1d = -gm * x / tembot;
        let y1d = -gm * y / tembot;
        let xd = x1;
        let yd = y1;

        x += h * xd;
        y += h * yd;
        x1 += h * x1d;
        y1 += h * y1d;
        t += h;

        // Second evaluation for RK2 correction
        let tembot = (x * x + y * y).powf(1.5);
        let x1d = -gm * x / tembot;
        let y1d = -gm * y / tembot;
        let xd = x1;
        let yd = y1;

        // RK2 averaging
        x = (xold + x) / 2.0 + 0.5 * h * xd;
        y = (yold + y) / 2.0 + 0.5 * h * yd;
        x1 = (x1old + x1) / 2.0 + 0.5 * h * x1d;
        y1 = (y1old + y1) / 2.0 + 0.5 * h * y1d;
    }

    PredictResult { xtf: x, ytf: y }
}

/// Lambert-guided missile trajectory prediction (predict45.m)
///
/// Predicts missile trajectory with Lambert guidance during boost phase.
///
/// # Arguments
/// * `tp` - Current time
/// * `xtp`, `ytp` - Current position
/// * `xtdp`, `ytdp` - Current velocity
/// * `tf` - Time to propagate to
/// * `tftot` - Total flight time
/// * `tupt` - Time for upward thrust
/// * `xf`, `yf` - Target position
/// * `itgt` - Target type (1 or 2, affects thrust profile)
///
/// # Returns
/// Final position (xtf, ytf)
#[allow(clippy::too_many_arguments)]
pub fn predict45(
    tp: f64,
    xtp: f64, ytp: f64,
    xtdp: f64, ytdp: f64,
    tf: f64,
    tftot: f64,
    tupt: f64,
    xf: f64, yf: f64,
    itgt: i32,
) -> PredictResult {
    let tpz = if itgt == 1 { 180.0 } else { 240.0 };
    let switch1 = 0;
    let gm = GM_FT;
    let h = 0.01;

    let mut t = tp;
    let mut xt = xtp;
    let mut yt = ytp;
    let zt = 0.0;
    let mut xtd = xtdp;
    let mut ytd = ytdp;
    let _ztd = 0.0;
    let zf = 0.0;

    let mut qboost = true;
    let mut axt = 0.0;
    let mut ayt = 0.0;

    while t <= tf - 0.00001 {
        let xtold = xt;
        let ytold = yt;
        let xtdold = xtd;
        let ytdold = ytd;

        // First Euler step
        let tembot = (xt * xt + yt * yt).powf(1.5);
        let mut xtdd = -gm * xt / tembot + axt;
        let mut ytdd = -gm * yt / tembot + ayt;

        // Compute thrust based on target type
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

        let atp = G_ACCEL * trst / wgt;

        xt += h * xtd;
        yt += h * ytd;
        xtd += h * xtdd;
        ytd += h * ytdd;
        t += h;

        // Second evaluation
        let tembot = (xt * xt + yt * yt).powf(1.5);
        xtdd = -gm * xt / tembot + axt;
        ytdd = -gm * yt / tembot + ayt;

        // RK2 averaging
        xt = (xtold + xt) / 2.0 + 0.5 * h * xtd;
        yt = (ytold + yt) / 2.0 + 0.5 * h * ytd;
        xtd = (xtdold + xtd) / 2.0 + 0.5 * h * xtdd;
        ytd = (ytdold + ytd) / 2.0 + 0.5 * h * ytdd;

        // Lambert guidance during boost
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
                qboost = true;
                axt = 0.0;
                ayt = 0.0;
                xtd = vtx;
                ytd = vty;
            } else {
                qboost = false;
                axt = 0.0;
                ayt = 0.0;
            }

            // Upward thrust phase
            if t < tupt {
                let rtmag = (xt * xt + yt * yt).sqrt();
                axt = atp * xt / rtmag;
                ayt = atp * yt / rtmag;
            }
        }
    }

    PredictResult { xtf: xt, ytf: yt }
}

/// Result from predictg with miss distance
#[derive(Debug, Clone, Copy)]
pub struct PredictGResult {
    pub xtf: f64,
    pub ytf: f64,
    pub zem1: f64,
    pub zem2: f64,
}

/// Gravity + thrust trajectory prediction with missile tracking (predictg.m)
///
/// Predicts trajectories of both interceptor and target, computing miss distance.
///
/// # Arguments
/// * `tdum` - Initial time
/// * `tf` - Final time
/// * `x`, `y` - Interceptor position
/// * `x1`, `y1` - Interceptor velocity
/// * `wp1` - Stage 1 propellant weight
/// * `wtot` - Total initial weight
/// * `tb1` - Stage 1 burn time
/// * `trst1` - Stage 1 thrust
/// * `tb2` - Stage 2 burn time
/// * `wp2` - Stage 2 propellant weight
/// * `wtot2` - Weight at stage 2 ignition
/// * `trst2` - Stage 2 thrust
/// * `wpay` - Payload weight
/// * `xm`, `ym` - Target position
/// * `x1m`, `y1m` - Target velocity
/// * `tgo` - Time to go
#[allow(clippy::too_many_arguments)]
pub fn predict_g(
    tdum: f64, tf: f64,
    x: f64, y: f64, x1: f64, y1: f64,
    wp1: f64, wtot: f64, tb1: f64, trst1: f64,
    tb2: f64, wp2: f64, wtot2: f64, trst2: f64,
    wpay: f64,
    xm: f64, ym: f64, x1m: f64, y1m: f64,
    tgo: f64,
) -> PredictGResult {
    let h = if tgo > 1.0 { 0.01 } else { tgo };
    let gm = GM_FT;

    let mut t = tdum;
    let mut x = x;
    let mut y = y;
    let mut x1 = x1;
    let mut y1 = y1;
    let mut xm = xm;
    let mut ym = ym;
    let mut x1m = x1m;
    let mut y1m = y1m;

    while t <= tf - 0.00001 {
        let xold = x;
        let yold = y;
        let x1old = x1;
        let y1old = y1;
        let xoldm = xm;
        let yoldm = ym;
        let x1oldm = x1m;
        let y1oldm = y1m;

        // Compute thrust and weight
        let (wgt, trst) = if t < tb1 {
            (-wp1 * t / tb1 + wtot, trst1)
        } else if t < tb1 + tb2 {
            (-wp2 * t / tb2 + wtot2 + wp2 * tb1 / tb2, trst2)
        } else {
            (wpay, 0.0)
        };

        let at = G_ACCEL * trst / wgt;
        let vel = (x1 * x1 + y1 * y1).sqrt();
        let axt = at * x1 / vel;
        let ayt = at * y1 / vel;

        // Interceptor dynamics
        let tembott = (x * x + y * y).powf(1.5);
        let x1d = -gm * x / tembott + axt;
        let y1d = -gm * y / tembott + ayt;
        let xd = x1;
        let yd = y1;

        // Target dynamics (ballistic)
        let tembotm = (xm * xm + ym * ym).powf(1.5);
        let x1dm = -gm * xm / tembotm;
        let y1dm = -gm * ym / tembotm;
        let xdm = x1m;
        let ydm = y1m;

        // Euler step
        x += h * xd;
        y += h * yd;
        x1 += h * x1d;
        y1 += h * y1d;
        xm += h * xdm;
        ym += h * ydm;
        x1m += h * x1dm;
        y1m += h * y1dm;
        t += h;

        // Second evaluation for RK2
        let vel = (x1 * x1 + y1 * y1).sqrt();
        let (wgt, trst) = if t < tb1 {
            (-wp1 * t / tb1 + wtot, trst1)
        } else if t < tb1 + tb2 {
            (-wp2 * t / tb2 + wtot2 + wp2 * tb1 / tb2, trst2)
        } else {
            (wpay, 0.0)
        };
        let at = G_ACCEL * trst / wgt;
        let axt = if vel > 0.0 { at * x1 / vel } else { 0.0 };
        let ayt = if vel > 0.0 { at * y1 / vel } else { 0.0 };

        let tembott = (x * x + y * y).powf(1.5);
        let x1d = -gm * x / tembott + axt;
        let y1d = -gm * y / tembott + ayt;
        let xd = x1;
        let yd = y1;

        let tembotm = (xm * xm + ym * ym).powf(1.5);
        let x1dm = -gm * xm / tembotm;
        let y1dm = -gm * ym / tembotm;
        let xdm = x1m;
        let ydm = y1m;

        // RK2 averaging
        x = (xold + x) / 2.0 + 0.5 * h * xd;
        y = (yold + y) / 2.0 + 0.5 * h * yd;
        x1 = (x1old + x1) / 2.0 + 0.5 * h * x1d;
        y1 = (y1old + y1) / 2.0 + 0.5 * h * y1d;
        xm = (xoldm + xm) / 2.0 + 0.5 * h * xdm;
        ym = (yoldm + ym) / 2.0 + 0.5 * h * ydm;
        x1m = (x1oldm + x1m) / 2.0 + 0.5 * h * x1dm;
        y1m = (y1oldm + y1m) / 2.0 + 0.5 * h * y1dm;
    }

    PredictGResult {
        xtf: x,
        ytf: y,
        zem1: x - xm,
        zem2: y - ym,
    }
}

/// Two-stage rocket trajectory prediction (predictb.m)
///
/// Propagates two-stage rocket position under gravity and thrust.
///
/// # Arguments
/// * `tf` - Time to propagate
/// * `x`, `y` - Initial position
/// * `x1`, `y1` - Initial velocity
/// * `wp1` - Stage 1 propellant weight
/// * `wtot` - Total initial weight
/// * `tb1` - Stage 1 burn time
/// * `trst1` - Stage 1 thrust
/// * `tb2` - Stage 2 burn time
/// * `wp2` - Stage 2 propellant weight
/// * `wtot2` - Weight at stage 2 ignition
/// * `trst2` - Stage 2 thrust
/// * `wpay` - Payload weight
#[allow(clippy::too_many_arguments)]
pub fn predict_b(
    tf: f64,
    xdum: f64, ydum: f64, x1dum: f64, y1dum: f64,
    wp1: f64, wtot: f64, tb1: f64, trst1: f64,
    tb2: f64, wp2: f64, wtot2: f64, trst2: f64,
    wpay: f64,
) -> PredictResult {
    let h = 0.01;
    let gm = GM_FT;

    let mut t = 0.0;
    let mut x = xdum;
    let mut y = ydum;
    let mut x1 = x1dum;
    let mut y1 = y1dum;

    while t <= tf - 0.00001 {
        let xold = x;
        let yold = y;
        let x1old = x1;
        let y1old = y1;

        // First Euler step
        let (wgt, trst) = if t < tb1 {
            (-wp1 * t / tb1 + wtot, trst1)
        } else if t < tb1 + tb2 {
            (-wp2 * t / tb2 + wtot2 + wp2 * tb1 / tb2, trst2)
        } else {
            (wpay, 0.0)
        };

        let at = G_ACCEL * trst / wgt;
        let vel = (x1 * x1 + y1 * y1).sqrt();
        let axt = at * x1 / vel;
        let ayt = at * y1 / vel;

        let tembott = (x * x + y * y).powf(1.5);
        let x1d = -gm * x / tembott + axt;
        let y1d = -gm * y / tembott + ayt;
        let xd = x1;
        let yd = y1;

        x += h * xd;
        y += h * yd;
        x1 += h * x1d;
        y1 += h * y1d;
        t += h;

        // Second evaluation for RK2
        let vel = (x1 * x1 + y1 * y1).sqrt();
        let (wgt, trst) = if t < tb1 {
            (-wp1 * t / tb1 + wtot, trst1)
        } else if t < tb1 + tb2 {
            (-wp2 * t / tb2 + wtot2 + wp2 * tb1 / tb2, trst2)
        } else {
            (wpay, 0.0)
        };
        let at = G_ACCEL * trst / wgt;
        let axt = if vel > 0.0 { at * x1 / vel } else { 0.0 };
        let ayt = if vel > 0.0 { at * y1 / vel } else { 0.0 };

        let tembott = (x * x + y * y).powf(1.5);
        let x1d = -gm * x / tembott + axt;
        let y1d = -gm * y / tembott + ayt;
        let xd = x1;
        let yd = y1;

        // RK2 averaging
        x = (xold + x) / 2.0 + 0.5 * h * xd;
        y = (yold + y) / 2.0 + 0.5 * h * yd;
        x1 = (x1old + x1) / 2.0 + 0.5 * h * x1d;
        y1 = (y1old + y1) / 2.0 + 0.5 * h * y1d;
    }

    PredictResult { xtf: x, ytf: y }
}

/// Calculate 3D great-circle distance on Earth (distance3dkm.m)
///
/// # Arguments
/// * `xt`, `yt`, `zt` - First position (feet)
/// * `xf`, `yf`, `zf` - Second position (feet)
///
/// # Returns
/// Distance in kilometers
pub fn distance3dkm(xt: f64, yt: f64, zt: f64, xf: f64, yf: f64, zf: f64) -> f64 {
    let r = (xt * xt + yt * yt + zt * zt).sqrt();
    let rf = (xf * xf + yf * yf + zf * zf).sqrt();
    let a = EARTH_RADIUS_FT;

    let cbeta = (xt * xf + yt * yf + zt * zf) / (r * rf);

    if cbeta <= 1.0 {
        let beta = cbeta.acos();
        a * beta / 3280.0
    } else {
        (xf - xt) / 3280.0
    }
}

/// Calculate 3D distance (distance3d.m) - returns in nmi
pub fn distance3d(xt: f64, yt: f64, zt: f64, xf: f64, yf: f64, zf: f64) -> f64 {
    let r = (xt * xt + yt * yt + zt * zt).sqrt();
    let rf = (xf * xf + yf * yf + zf * zf).sqrt();
    let a = EARTH_RADIUS_FT;

    let cbeta = (xt * xf + yt * yf + zt * zf) / (r * rf);

    if cbeta <= 1.0 {
        let beta = cbeta.acos();
        a * beta / 6076.0  // Convert to nautical miles
    } else {
        (xf - xt) / 6076.0
    }
}

/// Launch logic for 3D engagement (LAUNCHLOGIC.m)
///
/// Computes required missile velocity to intercept target
///
/// # Returns
/// (vm1, vm2, vm3, tf) - Required velocity components and flight time
#[allow(clippy::too_many_arguments)]
pub fn launch_logic(
    rm1: f64, rm2: f64, rm3: f64,
    rt1: f64, rt2: f64, rt3: f64,
    vt1: f64, vt2: f64, vt3: f64,
    vm: f64,
) -> (f64, f64, f64, f64) {
    let mut vm1 = 0.0;
    let mut vm2 = 0.0;
    let mut vm3 = 0.0;
    let mut tf_result = 0.1;

    let mut tf = 0.1;
    while tf <= 10.0 {
        let rt1f = rt1 + vt1 * tf;
        let rt2f = rt2 + vt2 * tf;
        let rt3f = rt3 + vt3 * tf;

        let arg = (rt2f - rm2) / (vm * tf);
        if arg.abs() <= 1.0 {
            let thet = arg.asin();
            let psi = (rt3f - rm3).atan2(rt1f - rm1);

            vm1 = vm * thet.cos() * psi.cos();
            vm2 = vm * thet.sin();
            vm3 = vm * thet.cos() * psi.sin();

            let rm1f = rm1 + vm1 * tf;
            let rm2f = rm2 + vm2 * tf;
            let rm3f = rm3 + vm3 * tf;

            let rtm1f = rt1f - rm1f;
            let rtm2f = rt2f - rm2f;
            let rtm3f = rt3f - rm3f;
            let rtmf = (rtm1f * rtm1f + rtm2f * rtm2f + rtm3f * rtm3f).sqrt();

            let vtm1 = vt1 - vm1;
            let vtm2 = vt2 - vm2;
            let vtm3 = vt3 - vm3;

            let vc = -(rtm1f * vtm1 + rtm2f * vtm2 + rtm3f * vtm3) / rtmf;

            if vc < 0.0 {
                tf_result = tf;
                break;
            }
        }
        tf_result = tf;
        tf += 0.1;
    }

    (vm1, vm2, vm3, tf_result)
}

/// Predict target position at future time (predict44.m)
///
/// Lambert-guided ballistic trajectory prediction using gravity-only propagation
/// after getting initial velocity from Lambert solver.
#[allow(clippy::too_many_arguments)]
pub fn predict44(
    tp: f64,
    xtp: f64, ytp: f64,
    _xtdp: f64, _ytdp: f64,
    tf: f64,
    tftot: f64,
    xf: f64, yf: f64,
) -> PredictResult {
    let switch1 = 0;
    let gm = GM_FT;
    let h = 0.01;

    let mut t = tp;
    let mut xt = xtp;
    let mut yt = ytp;
    let zt = 0.0;
    let zf = 0.0;

    // Get initial velocity from Lambert
    let tgolam = tftot - t;
    let result = lambert3d(xt, yt, zt, tgolam, xf, yf, zf, switch1);
    let mut xtd = result.vrx;
    let mut ytd = result.vry;

    while t <= tf - 0.00001 {
        let xtold = xt;
        let ytold = yt;
        let xtdold = xtd;
        let ytdold = ytd;

        // First derivative evaluation
        let tembot = (xt * xt + yt * yt).powf(1.5);
        let xtdd = -gm * xt / tembot;
        let ytdd = -gm * yt / tembot;

        // Euler step
        xt += h * xtd;
        yt += h * ytd;
        xtd += h * xtdd;
        ytd += h * ytdd;
        t += h;

        // Second evaluation
        let tembot = (xt * xt + yt * yt).powf(1.5);
        let xtdd = -gm * xt / tembot;
        let ytdd = -gm * yt / tembot;

        // RK2 averaging
        xt = (xtold + xt) / 2.0 + 0.5 * h * xtd;
        yt = (ytold + yt) / 2.0 + 0.5 * h * ytd;
        xtd = (xtdold + xtd) / 2.0 + 0.5 * h * xtdd;
        ytd = (ytdold + ytd) / 2.0 + 0.5 * h * ytdd;
    }

    PredictResult { xtf: xt, ytf: yt }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_predict_pz() {
        let a = EARTH_RADIUS_FT;
        let result = predict_pz(100.0, a, 0.0, 0.0, 25000.0);
        // Object should have moved
        assert!(result.ytf > 0.0);
    }

    #[test]
    fn test_distance3dkm() {
        let a = EARTH_RADIUS_FT;
        let dist = distance3dkm(a, 0.0, 0.0, a, 0.0, 0.0);
        assert!(dist.abs() < 0.001);
    }
}
