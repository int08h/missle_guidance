//! Chapter 44, Lesson 1: Ballistic Missile Trajectory
//!
//! Simulates a ballistic trajectory propagation using Lambert-derived initial velocity.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;
use crate::utils::{lambert3d, distance3dkm, EARTH_RADIUS_FT, GM_FT};

pub struct Results {
    pub time: Vec<f64>,
    pub dist_rt_knm: Vec<f64>,
    pub alt_t_km: Vec<f64>,
}

/// Run the C44L1 simulation
pub fn run() -> Results {
    let ts: f64 = 1.0;
    let rdeskm: f64 = 10000.0;
    let tloft: f64 = 0.0;
    let mut tftot = 252.0 + 0.223 * rdeskm - 5.44e-6 * rdeskm * rdeskm;
    tftot += tloft;

    let a = EARTH_RADIUS_FT;
    let gm = GM_FT;

    let xlongfdeg = 57.3 * rdeskm * 3280.0 / a;
    let xlongtdeg: f64 = 0.0;
    let xlatfdeg: f64 = 0.0;
    let xlattdeg: f64 = 0.0;
    let switch: i32 = 0;

    let mut t: f64 = 0.0;
    let mut s: f64 = 0.0;

    let xlongf = xlongfdeg / 57.3;
    let xlatf = xlatfdeg / 57.3;
    let xf = a * xlatf.cos() * xlongf.cos();
    let yf = a * xlatf.cos() * xlongf.sin();
    let zf: f64 = 0.0;

    let xlongt = xlongtdeg / 57.3;
    let xlatt = xlattdeg / 57.3;
    let mut xt = a * xlatt.cos() * xlongt.cos();
    let mut yt = a * xlatt.cos() * xlongt.sin();
    let zt: f64 = 0.0;
    let xtinit = xt;
    let ytinit = yt;
    let ztinit: f64 = 0.0;

    let h: f64 = 0.01;
    let mut alttkm = ((xt * xt + yt * yt).sqrt() - a) / 3280.0;

    // Get initial velocity from Lambert solver
    let tgolam = tftot - t;
    let result = lambert3d(xt, yt, zt, tgolam, xf, yf, zf, switch);
    let mut xtd = result.vrx;
    let mut ytd = result.vry;
    let _ztd = result.vrz;
    let _vbot = (xtd * xtd + ytd * ytd).sqrt() / 3280.0;

    let mut array_t = Vec::new();
    let mut array_distrtknm = Vec::new();
    let mut array_alttkm = Vec::new();

    while alttkm > -1.0 {
        let xtold = xt;
        let ytold = yt;
        let xtdold = xtd;
        let ytdold = ytd;

        // First derivative evaluation
        let tempbott = (xt * xt + yt * yt).powf(1.5);
        let xtdd = -gm * xt / tempbott;
        let ytdd = -gm * yt / tempbott;
        alttkm = ((xt * xt + yt * yt).sqrt() - a) / 3280.0;

        // Euler step
        xt += h * xtd;
        yt += h * ytd;
        xtd += h * xtdd;
        ytd += h * ytdd;
        t += h;

        // Second derivative evaluation (for RK2)
        let tempbott2 = (xt * xt + yt * yt).powf(1.5);
        let xtdd2 = -gm * xt / tempbott2;
        let ytdd2 = -gm * yt / tempbott2;

        // RK2 averaging
        xt = 0.5 * (xtold + xt + h * xtd);
        yt = 0.5 * (ytold + yt + h * ytd);
        xtd = 0.5 * (xtdold + xtd + h * xtdd2);
        ytd = 0.5 * (ytdold + ytd + h * ytdd2);

        s += h;

        if s >= ts - 0.0001 {
            s = 0.0;
            alttkm = ((xt * xt + yt * yt).sqrt() - a) / 3280.0;
            let distrtknm = distance3dkm(xt, yt, zt, xtinit, ytinit, ztinit);

            array_t.push(t);
            array_distrtknm.push(distrtknm);
            array_alttkm.push(alttkm);
        }
    }

    Results {
        time: array_t,
        dist_rt_knm: array_distrtknm,
        alt_t_km: array_alttkm,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c44l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.dist_rt_knm.clone(),
        results.alt_t_km.clone(),
    ])?;

    let plot_file = format!("{}/c44l1_trajectory.png", output_dir);
    let config = PlotConfig::new("Ballistic Missile Trajectory")
        .with_labels("Downrange (km)", "Altitude (km)");

    let series = vec![
        Series::new(results.dist_rt_knm.clone(), results.alt_t_km.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C44L1: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c44l1_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
