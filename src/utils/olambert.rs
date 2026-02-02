//! olambert.m port - Orbital Lambert solver with iteration history
//!
//! Solves Lambert's problem using a for-loop search approach and records
//! all iteration values, matching the MATLAB implementation exactly.

use std::f64::consts::PI;

/// Result from the olambert solver containing iteration history
#[derive(Debug, Clone)]
pub struct OLambertResult {
    pub vrx: f64,
    pub vry: f64,
    pub array_vrx: Vec<f64>,
    pub array_vry: Vec<f64>,
    pub array_tf: Vec<f64>,
    pub count: usize,
}

/// Solves the 2D Lambert problem using for-loop search (matching MATLAB olambert.m)
///
/// # Arguments
/// * `xic`, `yic` - Initial position
/// * `tfdes` - Desired time of flight
/// * `xf`, `yf` - Final position
/// * `xlongm` - Initial longitude (radians)
/// * `xlongt` - Target longitude (radians)
/// * `a` - Earth radius (ft)
/// * `gm` - Gravitational parameter (ft^3/s^2)
///
/// # Returns
/// `OLambertResult` containing final velocity and iteration history arrays
pub fn olambert(
    xic: f64, yic: f64,
    tfdes: f64,
    xf: f64, yf: f64,
    xlongm: f64, xlongt: f64,
    _a: f64,
    gm: f64,
) -> OLambertResult {
    // Initialize outputs
    let mut count: usize = 0;
    let mut vrx: f64 = 0.0;
    let mut vry: f64 = 0.0;
    let mut tf: f64;

    let mut array_vrx = Vec::new();
    let mut array_vry = Vec::new();
    let mut array_tf = Vec::new();

    // Section A
    let ric = (xic * xic + yic * yic).sqrt();
    let rf = (xf * xf + yf * yf).sqrt();
    let cphi = (xic * xf + yic * yf) / (ric * rf);
    let phi = cphi.acos();
    let r0 = ric;
    let degrad = 360.0 / (2.0 * PI);

    // Initialize while loop
    let mut second_time_through = 0;
    let mut gamdegnew: f64 = 0.0;
    let mut gamdegfin: f64 = 0.0;
    let mut gamdeg_break: f64 = 0.0;

    // Program executes this loop twice
    while second_time_through <= 1 {
        // Initialize for loop
        let (start, step, stop) = if second_time_through == 0 {
            (-90.0, 0.1, 90.0)
        } else {
            (gamdegnew, 0.0001, gamdegfin)
        };

        // Main body of program - for loop
        let mut gamdeg = start;
        while gamdeg <= stop {
            // Section B
            let gam = gamdeg / degrad;
            let top = gm * (1.0 - phi.cos());
            let temp = r0 * gam.cos() / rf - (phi + gam).cos();
            let bot = r0 * gam.cos() * temp;

            // Condition #1: not (top < 0 or bot < 0)
            if !(top < 0.0 || bot < 0.0) {
                // Section C
                let v = (top / bot).sqrt();
                if xlongt > xlongm {
                    vrx = v * (PI / 2.0 - gam + xlongm).cos();
                    vry = v * (PI / 2.0 - gam + xlongm).sin();
                } else {
                    vrx = v * (-PI / 2.0 + gam + xlongm).cos();
                    vry = v * (-PI / 2.0 + gam + xlongm).sin();
                }

                let xlam = r0 * v * v / gm;
                let top1 = gam.tan() * (1.0 - phi.cos()) + (1.0 - xlam) * phi.sin();
                let bot1p = (1.0 - phi.cos()) / (xlam * gam.cos() * gam.cos());
                let bot1 = (2.0 - xlam) * (bot1p + (gam + phi).cos() / gam.cos());
                let top2 = 2.0 * gam.cos();

                // Condition #2: not ((2/xlam - 1) < 0)
                if !((2.0 / xlam - 1.0) < 0.0) {
                    // Section D
                    let bot2 = xlam * (2.0 / xlam - 1.0).powf(1.5);
                    let top3 = (2.0 / xlam - 1.0).sqrt();
                    let bot3 = gam.cos() / (phi / 2.0).tan() - gam.sin();
                    let temp_val = (top2 / bot2) * top3.atan2(bot3);
                    tf = r0 * (top1 / bot1 + temp_val) / (v * gam.cos());

                    // Condition #3: break if tf > tfdes
                    if tf > tfdes {
                        gamdeg_break = gamdeg;
                        break;
                    }

                    // Output arrays (record iteration)
                    count += 1;
                    array_vrx.push(vrx);
                    array_vry.push(vry);
                    array_tf.push(tf);
                }
            }

            gamdeg += step;
            gamdeg_break = gamdeg;
        }

        // Section E
        gamdegnew = gamdeg_break - 0.15;
        gamdegfin = gamdeg_break + 1.0;

        second_time_through += 1;
    }

    OLambertResult {
        vrx,
        vry,
        array_vrx,
        array_vry,
        array_tf,
        count,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_olambert_basic() {
        let a = 2.0926e7;
        let gm = 1.4077e16;
        let degrad = 360.0 / (2.0 * PI);

        let xlongmdeg = 45.0;
        let xlongtdeg = 90.0;
        let xlongm = xlongmdeg / degrad;
        let xlongt = xlongtdeg / degrad;

        let xm = a * xlongm.cos();
        let ym = a * xlongm.sin();
        let xt = a * xlongt.cos();
        let yt = a * xlongt.sin();

        let result = olambert(xm, ym, 1000.0, xt, yt, xlongm, xlongt, a, gm);

        assert!(result.count > 0);
        assert!(!result.array_vrx.is_empty());
    }
}
