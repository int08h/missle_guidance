//! 2D Lambert trajectory solver (lamberpz.m)
//!
//! Solves Lambert's problem in 2D to find the initial velocity vector required to transfer
//! from an initial position to a final position in a specified time.

use std::f64::consts::PI;
use super::constants::{EARTH_RADIUS_FT, GM_FT};

/// Result from the Lambert 2D solver containing velocity components
#[derive(Debug, Clone, Copy)]
pub struct Lambert2DResult {
    pub vrx: f64,
    pub vry: f64,
}

/// Solves the 2D Lambert problem
///
/// # Arguments
/// * `xic`, `yic` - Initial position
/// * `tfdes` - Desired time of flight
/// * `xf`, `yf` - Final position
/// * `xlongm` - Initial longitude (radians)
/// * `xlongt` - Target longitude (radians)
///
/// # Returns
/// `Lambert2DResult` containing required velocity components (vrx, vry)
pub fn lambert2d(
    xic: f64, yic: f64,
    tfdes: f64,
    xf: f64, yf: f64,
    xlongm: f64, xlongt: f64,
) -> Lambert2DResult {
    let _a = EARTH_RADIUS_FT;
    let gm = GM_FT;

    let ric = (xic * xic + yic * yic).sqrt();
    let rf = (xf * xf + yf * yf).sqrt();
    let cphi = (xic * xf + yic * yf) / (ric * rf);
    let phi = cphi.acos();
    let sphi = phi.sin();
    let r0 = ric;

    let mut icount = 0;

    // Initial bounds for gamma angle
    let gmin_init = ((sphi - (2.0 * r0 * (1.0 - cphi) / rf).sqrt()) / (1.0 - cphi)).atan2(1.0);
    let gmax_init = ((sphi + (2.0 * r0 * (1.0 - cphi) / rf).sqrt()) / (1.0 - cphi)).atan2(1.0);

    let mut gmin = gmin_init;
    let mut gmax = gmax_init;
    let mut gam = (gmin + gmax) / 2.0;
    let mut tf = 0.0;
    let mut gold = gam;
    let mut told = 0.0;
    let mut vrx = 0.0;
    let mut vry = 0.0;

    while (tfdes - tf).abs() > 0.00000001 * tfdes {
        let top = gm * (1.0 - cphi);
        let temp = r0 * gam.cos() / rf - (phi + gam).cos();
        let bot = r0 * gam.cos() * temp;

        if bot <= 0.0 {
            // Invalid configuration, adjust bounds
            gmax = gam;
            gam = (gmin + gmax) / 2.0;
            continue;
        }

        let v = (top / bot).sqrt();

        // Compute velocity based on direction
        if xlongt > xlongm {
            vrx = v * (PI / 2.0 - gam + xlongm).cos();
            vry = v * (PI / 2.0 - gam + xlongm).sin();
        } else {
            vrx = v * (-PI / 2.0 + gam + xlongm).cos();
            vry = v * (-PI / 2.0 + gam + xlongm).sin();
        }

        let xlam = r0 * v * v / gm;
        let top1 = gam.tan() * (1.0 - cphi) + (1.0 - xlam) * sphi;
        let bot1p = (1.0 - cphi) / (xlam * gam.cos() * gam.cos());
        let bot1 = (2.0 - xlam) * (bot1p + (gam + phi).cos() / gam.cos());
        let top2 = 2.0 * gam.cos();

        let inner = 2.0 / xlam - 1.0;
        if inner < 0.0 {
            // Hyperbolic case - adjust bounds
            gmax = gam;
            gam = (gmin + gmax) / 2.0;
            continue;
        }

        let bot2 = xlam * inner.powf(1.5);
        let top3 = inner.sqrt();
        let bot3 = gam.cos() / (phi / 2.0).tan() - gam.sin();
        let temp_tf = (top2 / bot2) * top3.atan2(bot3);
        tf = r0 * (top1 / bot1 + temp_tf) / (v * gam.cos());

        icount += 1;

        if tf > tfdes {
            gmax = gam;
        } else {
            gmin = gam;
        }

        let xnext = if icount == 1 {
            (gmax + gmin) / 2.0
        } else {
            let mut xnext = gam + (gam - gold) * (tfdes - tf) / (tf - told);
            if xnext > gmax || xnext < gmin {
                xnext = (gmax + gmin) / 2.0;
            }
            xnext
        };

        gold = gam;
        told = tf;
        gam = xnext;

        if icount > 100 {
            break;
        }
    }

    Lambert2DResult { vrx, vry }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lambert2d_basic() {
        let a = EARTH_RADIUS_FT;
        let xic = a;
        let yic = 0.0;
        let xf = a * 0.8;
        let yf = a * 0.6;
        let tfdes = 1000.0;

        let result = lambert2d(xic, yic, tfdes, xf, yf, 0.0, 0.5);

        let v_mag = (result.vrx.powi(2) + result.vry.powi(2)).sqrt();
        assert!(v_mag > 0.0);
    }
}
