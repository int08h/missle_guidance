//! Chapter 40, Lesson 2: Lambert Trajectory Simulation
//!
//! Ballistic trajectory with Lambert guidance to specified endpoint.
//!
//! MATLAB NOTE: Calls `LAMBERT3D` and `distance3d` (which should be `distance3dkm`).
//! The function naming is inconsistent but the implementations are available.
//! Output: 3 columns [T, DISTRTNM, ALTNM]

use crate::save_data;
use crate::utils::{lambert3d, distance3d, EARTH_RADIUS_FT, GM_FT};

pub struct Results {
    pub time: Vec<f64>,
    pub dist_rt_nm: Vec<f64>,
    pub alt_nm: Vec<f64>,
}

/// Run the C40L2 simulation
pub fn run() -> Results {
    let switch1 = 0;
    let xlongtdeg: f64 = 7.42;
    let xlattdeg: f64 = 43.75;
    let xlatfdeg: f64 = 36.175;
    let xlongfdeg: f64 = -115.136;
    let tf: f64 = 2000.0;
    let a = EARTH_RADIUS_FT;
    let gm = GM_FT;
    let w: f64 = -6.283185 / 86400.0;
    let h: f64 = 0.001;

    let xlongf = xlongfdeg / 57.3 - w * tf;
    let xlatf = xlatfdeg / 57.3;
    let xlongt = xlongtdeg / 57.3;
    let xlatt = xlattdeg / 57.3;

    let xf = a * xlatf.cos() * xlongf.cos();
    let yf = a * xlatf.cos() * xlongf.sin();
    let zf = a * xlatf.sin();

    let mut xt = a * xlatt.cos() * xlongt.cos();
    let mut yt = a * xlatt.cos() * xlongt.sin();
    let mut zt = a * xlatt.sin();

    let xtinit = xt;
    let ytinit = yt;
    let ztinit = zt;

    // Get initial velocity from Lambert solver
    let result = lambert3d(xt, yt, zt, tf, xf, yf, zf, switch1);
    let mut xtd = result.vrx;
    let mut ytd = result.vry;
    let mut ztd = result.vrz;

    let mut t = 0.0;
    let mut s = 0.0;
    let mut altnm = ((xt * xt + yt * yt + zt * zt).sqrt() - a) / 6076.0;

    let mut array_t = Vec::new();
    let mut array_distrtnm = Vec::new();
    let mut array_altnm = Vec::new();

    while altnm > -1.0 {
        let xtold = xt;
        let ytold = yt;
        let ztold = zt;
        let xtdold = xtd;
        let ytdold = ytd;
        let ztdold = ztd;

        // First derivative evaluation
        let tempbott = (xt * xt + yt * yt + zt * zt).powf(1.5);
        let xtdd = -gm * xt / tempbott;
        let ytdd = -gm * yt / tempbott;
        let ztdd = -gm * zt / tempbott;

        // Euler step
        xt += h * xtd;
        yt += h * ytd;
        zt += h * ztd;
        xtd += h * xtdd;
        ytd += h * ytdd;
        ztd += h * ztdd;
        t += h;

        // Second evaluation
        let tempbott = (xt * xt + yt * yt + zt * zt).powf(1.5);
        let xtdd = -gm * xt / tempbott;
        let ytdd = -gm * yt / tempbott;
        let ztdd = -gm * zt / tempbott;
        altnm = ((xt * xt + yt * yt + zt * zt).sqrt() - a) / 6076.0;

        // RK2 averaging
        xt = 0.5 * (xtold + xt + h * xtd);
        yt = 0.5 * (ytold + yt + h * ytd);
        zt = 0.5 * (ztold + zt + h * ztd);
        xtd = 0.5 * (xtdold + xtd + h * xtdd);
        ytd = 0.5 * (ytdold + ytd + h * ytdd);
        ztd = 0.5 * (ztdold + ztd + h * ztdd);

        s += h;
        if s >= 9.9999 {
            s = 0.0;
            let xte = xt * (w * t).cos() - yt * (w * t).sin();
            let yte = xt * (w * t).sin() + yt * (w * t).cos();
            let zte = zt;

            let distrtnm = distance3d(xte, yte, zte, xtinit, ytinit, ztinit);
            altnm = ((xt * xt + yt * yt + zt * zt).sqrt() - a) / 6076.0;

            array_t.push(t);
            array_distrtnm.push(distrtnm);
            array_altnm.push(altnm);
        }
    }

    Results {
        time: array_t,
        dist_rt_nm: array_distrtnm,
        alt_nm: array_altnm,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c40l2_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.dist_rt_nm.clone(),
        results.alt_nm.clone(),
    ])?;

    println!("C40L2: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c40l2_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
