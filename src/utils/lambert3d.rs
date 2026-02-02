//! 3D Lambert trajectory solver
//!
//! Solves Lambert's problem to find the initial velocity vector required to transfer
//! from an initial position to a final position in a specified time.

use std::f64::consts::PI;
use super::constants::{GM_FT, HALF_PI};

#[cfg(test)]
use super::constants::EARTH_RADIUS_FT;

/// Result from the Lambert 3D solver containing velocity components
#[derive(Debug, Clone, Copy)]
pub struct Lambert3DResult {
    pub vrx: f64,
    pub vry: f64,
    pub vrz: f64,
}

/// Solves the 3D Lambert problem
///
/// # Arguments
/// * `xt`, `yt`, `zt` - Initial position (from position)
/// * `tf` - Time of flight
/// * `xf`, `yf`, `zf` - Final position (to position)
/// * `switch1` - 0 for short path, 1 for long path
///
/// # Returns
/// `Lambert3DResult` containing required velocity components (vrx, vry, vrz)
#[allow(clippy::too_many_arguments)]
pub fn lambert3d(
    xt: f64, yt: f64, zt: f64,
    tf: f64,
    xf: f64, yf: f64, zf: f64,
    switch1: i32,
) -> Lambert3DResult {
    let gm = GM_FT;

    // Relative position vector
    let rf0x = xf - xt;
    let rf0y = yf - yt;
    let rf0z = zf - zt;

    // Dot products
    let r0_dot_rf = xt * xf + yt * yf + zt * zf;
    let r0_dot_rf0 = xt * rf0x + yt * rf0y + zt * rf0z;

    // Magnitudes
    let r0_mag = (xt * xt + yt * yt + zt * zt).sqrt();
    let rf_mag = (xf * xf + yf * yf + zf * zf).sqrt();
    let rf0_mag = (rf0x * rf0x + rf0y * rf0y + rf0z * rf0z).sqrt();

    let ratio = r0_mag / rf_mag;
    let gm_div_r0 = gm / r0_mag;

    let cos_t = r0_dot_rf / (r0_mag * rf_mag);
    let vnumer = gm_div_r0 * (1.0 - cos_t);

    // Set bounds for gamma angle based on path selection
    let (mut g_min, mut g_max, theta) = if switch1 == 0 {
        let g_min = HALF_PI - (r0_dot_rf0 / (r0_mag * rf0_mag)).acos();
        let g_max = HALF_PI;
        let theta = cos_t.acos();
        (g_min, g_max, theta)
    } else {
        let g_min = -HALF_PI;
        let g_max = -HALF_PI + (r0_dot_rf0 / (r0_mag * rf0_mag)).acos();
        let theta = 2.0 * PI - cos_t.acos();
        (g_min, g_max, theta)
    };

    let sin_t = theta.sin();
    let cot_half_t = 1.0 / (theta / 2.0).tan();

    let mut gamma = (g_max + g_min) / 2.0;
    let mut gold = g_min;
    let mut told = 0.0;
    let mut t = 0.0;
    let mut t_min = 0.0;
    let mut iters = 1;
    let mut rv0_mag = 0.0;

    // Semi-perimeter for normalization
    let s = 0.5 * (r0_mag + rf_mag + rf0_mag);
    let _bl = (r0_mag * rf_mag).sqrt() * (theta / 2.0).cos() / s;
    let _bt = (8.0 * gm / (s * s * s)).sqrt() * tf;

    // Iterative solution using binary search with secant acceleration
    while (tf - t).abs() > 0.00000001 * tf {
        let sin_g = gamma.sin();
        let cos_g = gamma.cos();
        let tan_g = sin_g / cos_g;
        let cos_t_plus_g = (theta + gamma - 2.0 * PI).cos();
        let term1 = (ratio * cos_g - cos_t_plus_g) * cos_g;
        rv0_mag = (vnumer / term1).sqrt();
        let lambda = rv0_mag * rv0_mag / gm_div_r0;

        if lambda < 1.9999999 {
            // Elliptic case
            let term0 = (2.0 / lambda - 1.0).sqrt();
            let term1_t = (tan_g * (1.0 - cos_t) + (1.0 - lambda) * sin_t)
                / ((2.0 - lambda) * ratio);
            let term2 = (cos_g + cos_g) / (lambda * term0 * term0 * term0);
            let term3 = (term0).atan2(cos_g * cot_half_t - sin_g);
            t = (r0_mag / (rv0_mag * cos_g)) * (term1_t + term2 * term3);
        } else if lambda > 2.0000001 {
            // Hyperbolic case
            let term0 = (1.0 - 2.0 / lambda).sqrt();
            let term1_t = (tan_g * (1.0 - cos_t) + (1.0 - lambda) * sin_t)
                / ((2.0 - lambda) * ratio);
            let term2 = cos_g / (lambda * term0 * term0 * term0);
            let term3_base = sin_g - cos_g * cot_half_t;
            let term3 = ((term3_base - term0) / (term3_base + term0)).ln();
            t = (r0_mag / (rv0_mag * cos_g)) * (term1_t - term2 * term3);
        } else {
            // Parabolic case
            let term0 = cos_g * cot_half_t;
            let term1_p = term0 - sin_g;
            let term0_p = (3.0 * term0 * term1_p + 1.0) / (term1_p * term1_p * term1_p);
            t = term0_p * (2.0 * r0_mag) / (3.0 * rv0_mag);
        }

        // Update bounds
        if t > tf && gamma < g_max {
            g_max = gamma;
        }
        if t < 0.0 && gamma < g_max {
            g_max = gamma;
        }
        if t < tf && gamma > g_min {
            g_min = gamma;
            t_min = t;
        }

        // Compute next gamma using secant method or bisection
        let next = if t < 0.0 {
            gold = g_min;
            told = t_min;
            (g_min + g_max) / 2.0
        } else {
            let mut next = gamma + (tf - t) * (gamma - gold) / (t - told);
            if next >= g_max {
                next = (gamma + g_max) / 2.0;
            } else if next <= g_min {
                next = (gamma + g_min) / 2.0;
            }
            gold = gamma;
            told = t;
            next
        };

        gamma = next;
        iters += 1;
        if iters > 100 {
            break;
        }
    }

    // Compute velocity direction
    let (vunit_x, vunit_y, vunit_z) = if switch1 == 0 {
        let angle = HALF_PI - gamma;
        let sina = angle.sin();
        let cosa = angle.cos();

        let v1x = xt;
        let v1y = yt;
        let v1z = zt;
        let v2x = xf;
        let v2y = yf;
        let v2z = zf;

        let mag1 = (v1x * v1x + v1y * v1y + v1z * v1z).sqrt();
        let dotmag = v1x * v2x + v1y * v2y + v1z * v2z;
        let crossx = v1y * v2z - v1z * v2y;
        let crossy = v1z * v2x - v1x * v2z;
        let crossz = v1x * v2y - v1y * v2x;
        let crossmag = (crossx * crossx + crossy * crossy + crossz * crossz).sqrt();

        let c2 = mag1 * sina / crossmag;
        let c1 = cosa / mag1 - dotmag * c2 / (mag1 * mag1);

        let rtempx = c1 * v1x;
        let rtempy = c1 * v1y;
        let rtempz = c1 * v1z;

        let vunitx = c2 * v2x + rtempx;
        let vunity = c2 * v2y + rtempy;
        let vunitz = c2 * v2z + rtempz;

        (vunitx, vunity, vunitz)
    } else {
        let angle = gamma - HALF_PI;
        let sina = angle.sin();
        let cosa = angle.cos();

        let v1x = xt;
        let v1y = yt;
        let v1z = zt;
        let v2x = xf;
        let v2y = yf;
        let v2z = zf;

        let mag1 = (v1x * v1x + v1y * v1y + v1z * v1z).sqrt();
        let dotmag = v1x * v2x + v1y * v2y + v1z * v2z;
        let crossx = v1y * v2z - v1z * v2y;
        let crossy = v1z * v2x - v1x * v2z;
        let crossz = v1x * v2y - v1y * v2x;
        let crossmag = (crossx * crossx + crossy * crossy + crossz * crossz).sqrt();

        let c2 = mag1 * sina / crossmag;
        let c1 = cosa / mag1 - dotmag * c2 / (mag1 * mag1);

        let rtempx = c1 * v1x;
        let rtempy = c1 * v1y;
        let rtempz = c1 * v1z;

        let vunitx = c2 * v2x + rtempx;
        let vunity = c2 * v2y + rtempy;
        let vunitz = c2 * v2z + rtempz;

        (vunitx, vunity, vunitz)
    };

    Lambert3DResult {
        vrx: rv0_mag * vunit_x,
        vry: rv0_mag * vunit_y,
        vrz: rv0_mag * vunit_z,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lambert3d_basic() {
        // Test with typical orbital mechanics values
        let a = EARTH_RADIUS_FT;
        let xt = a;
        let yt = 0.0;
        let zt = 0.0;
        let xf = 0.0;
        let yf = a * 1.5;
        let zf = 0.0;
        let tf = 1000.0;

        let result = lambert3d(xt, yt, zt, tf, xf, yf, zf, 0);

        // Velocity should be non-zero and reasonable
        let v_mag = (result.vrx.powi(2) + result.vry.powi(2) + result.vrz.powi(2)).sqrt();
        assert!(v_mag > 0.0);
        assert!(v_mag < 100000.0); // Reasonable orbital velocity
    }
}
