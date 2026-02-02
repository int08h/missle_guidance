//! Chapter 43, Lesson 1: Ballistic Missile Trajectory with Lambert Guidance
//!
//! Simulates a ballistic missile trajectory using Lambert guidance during boost phase.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;
use crate::utils::{lambert3d, distance3dkm, EARTH_RADIUS_FT, GM_FT};

pub struct Results {
    pub time: Vec<f64>,
    pub dist_nm: Vec<f64>,
    pub alt_nm: Vec<f64>,
    pub q: Vec<f64>,
    pub isee: Vec<f64>,
}

/// Run the C43L1 simulation
pub fn run() -> Results {
    let rdeskm: f64 = 7000.0;
    #[allow(unused_assignments)]
    let mut tf: f64 = 0.0;
    let tfinish: f64 = 999999.0;
    let tloft: f64 = 500.0;
    let tgravend: f64 = 100.0;
    let gamdegic: f64 = 89.8;
    let tupt: f64 = 20.0;
    let rdesrkm: f64 = 560.0;
    let switch: i32 = 0;

    let mut qboost: bool = true;
    let mut qfirst: bool = true;
    let gamdeg = gamdegic;

    let h: f64 = 0.01;
    let mut t: f64 = 0.0;
    let mut s: f64 = 0.0;

    let a = EARTH_RADIUS_FT;
    let gm = GM_FT;

    let altnm: f64 = 0.0;
    let mut alt = altnm * 6076.0;
    let angdeg: f64 = 0.0;
    let ang = angdeg / 57.3;

    let mut x = (a + alt) * ang.cos();
    let mut y = (a + alt) * ang.sin();
    // alt computed fresh in the loop before use
    let xfirst = x;
    let yfirst = y;
    let zfirst: f64 = 0.0;

    let mut x1 = (1.5708 - gamdeg / 57.3 + ang).cos();
    let mut y1 = (1.5708 - gamdeg / 57.3 + ang).sin();
    let mut axt: f64 = 0.0;
    let mut ayt: f64 = 0.0;

    let xlongtdeg = 57.3 * rdeskm * 3280.0 / a;
    let xlongrdeg = 57.3 * rdesrkm * 3280.0 / a;
    tf = 252.0 + 0.223 * rdeskm - 5.44e-6 * rdeskm * rdeskm;
    tf += tloft;
    let xlongt = xlongtdeg / 57.3;
    let xlongr = xlongrdeg / 57.3;
    let xf = a * xlongt.cos();
    let yf = a * xlongt.sin();
    let xr = a * xlongr.cos();
    let yr = a * xlongr.sin();

    let z: f64 = 0.0;
    let zf: f64 = 0.0;

    let mut array_t = Vec::new();
    let mut array_distnm = Vec::new();
    let mut array_altnm = Vec::new();
    let mut array_q = Vec::new();
    let mut array_isee = Vec::new();

    let mut altnm_val = ((x * x + y * y).sqrt() - a) / 3280.0;

    while !(altnm_val < -1.0 || t > tfinish) {
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
        x += h * xd;
        y += h * yd;
        x1 += h * x1d;
        y1 += h * y1d;
        t += h;

        // Second derivative evaluation (for RK2)
        let xd2 = x1;
        let yd2 = y1;

        // RK2 averaging
        x = (xold + x) / 2.0 + 0.5 * h * xd2;
        y = (yold + y) / 2.0 + 0.5 * h * yd2;
        x1 = (x1old + x1) / 2.0 + 0.5 * h * x1d;
        y1 = (y1old + y1) / 2.0 + 0.5 * h * y1d;

        s += h;
        let tgolam = tf - t;

        // Lambert guidance during boost
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

        // Sampling
        if s >= 0.9999 {
            s = 0.0;
            let distnm = distance3dkm(x, y, z, xfirst, yfirst, zfirst);
            altnm_val = ((x * x + y * y).sqrt() - a) / 3280.0;
            let rmag = (x * x + y * y).sqrt();
            let vmag = (xd * xd + yd * yd).sqrt();
            let _gamdeg = 90.0 - 57.3 * ((x * xd + y * yd) / (rmag * vmag)).acos();
            let rho = 0.0034 * (-alt / 22000.0).exp();
            let q_val = 0.5 * rho * vel * vel;
            let rrmag = (xr * xr + yr * yr).sqrt();
            let rrtmag = ((x - xr).powi(2) + (y - yr).powi(2)).sqrt();
            let eldeg = 90.0 - 57.3 * ((xr * (x - xr) + yr * (y - yr)) / (rrmag * rrtmag)).acos();

            let isee = if eldeg > 2.0 && eldeg < 85.0 { 1.0 } else { 0.0 };

            array_t.push(t);
            array_distnm.push(distnm);
            array_altnm.push(altnm_val);
            array_q.push(q_val);
            array_isee.push(isee);
        }
    }

    Results {
        time: array_t,
        dist_nm: array_distnm,
        alt_nm: array_altnm,
        q: array_q,
        isee: array_isee,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c43l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.dist_nm.clone(),
        results.alt_nm.clone(),
        results.q.clone(),
        results.isee.clone(),
    ])?;

    let plot_file = format!("{}/c43l1_trajectory.png", output_dir);
    let config = PlotConfig::new("Ballistic Missile Trajectory")
        .with_labels("Downrange (km)", "Altitude (km)");

    let series = vec![
        Series::new(results.dist_nm.clone(), results.alt_nm.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C43L1: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c43l1_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
