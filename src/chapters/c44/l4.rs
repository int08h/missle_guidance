//! Chapter 44, Lesson 4: BMD Target Coverage Analysis
//!
//! Ballistic missile defense: coverage analysis for varying target locations.

use crate::save_data;
use crate::utils::{lambert3d, kepler1, StateVector, distance3dkm, EARTH_RADIUS_FT, GM_FT};

pub struct Results {
    pub tf: Vec<f64>,
    pub tlaunch: Vec<f64>,
    pub xlongt: Vec<f64>,
    pub xlatt: Vec<f64>,
    pub vbom: Vec<f64>,
    pub vbot: Vec<f64>,
}

/// Run the C44L4 simulation
pub fn run() -> Results {
    let tlaunch_start: f64 = 300.0;
    let rdeskm: f64 = 10000.0;
    let altmkmic: f64 = 0.0;
    let tloft: f64 = 0.0;
    let vbolim: f64 = 5.0;
    let xlongmdeg: f64 = 60.0;
    let xlatmdeg: f64 = 0.0;
    let switch: i32 = 0;
    let switchm: i32 = 0;

    let a = EARTH_RADIUS_FT;
    let _gm = GM_FT;

    let xlongfdeg = 57.3 * rdeskm * 3280.0 / a;
    let xlatfdeg: f64 = 0.0;

    let xlongf = xlongfdeg / 57.3;
    let xlatf = xlatfdeg / 57.3;
    let xlongm = xlongmdeg / 57.3;
    let xlatm = xlatmdeg / 57.3;

    let xf = a * xlatf.cos() * xlongf.cos();
    let yf = a * xlatf.cos() * xlongf.sin();
    let zf = a * xlatf.sin();

    let altm = altmkmic * 3280.0;
    let xm = (a + altm) * xlatm.cos() * xlongm.cos();
    let ym = (a + altm) * xlatm.cos() * xlongm.sin();
    let zm = (a + altm) * xlatm.sin();

    let mut array_tf = Vec::new();
    let mut array_tlaunch = Vec::new();
    let mut array_xlongt = Vec::new();
    let mut array_xlatt = Vec::new();
    let mut array_vbom = Vec::new();
    let mut array_vbot = Vec::new();

    let mut xlongtdeg: f64 = -100.0;
    while xlongtdeg <= 200.0 {
        let mut xlattdeg: f64 = -60.0;
        while xlattdeg <= 60.0 {
            let mut tlaunch = tlaunch_start;
            while tlaunch <= 1600.0 {
                let xlongt = xlongtdeg / 57.3;
                let xlatt = xlattdeg / 57.3;

                let xt = a * xlatt.cos() * xlongt.cos();
                let yt = a * xlatt.cos() * xlongt.sin();
                let zt = a * xlatt.sin();

                let distfkm = distance3dkm(xf, yf, zf, xt, yt, zt);
                let mut tftot = 252.0 + 0.223 * distfkm - 5.44e-6 * distfkm * distfkm;
                tftot += tloft;

                let mut tf = tlaunch + 60.0;
                while tf <= tftot {
                    let tgolam = tftot;
                    let result = lambert3d(xt, yt, zt, tgolam, xf, yf, zf, switch);
                    let xtd = result.vrx;
                    let ytd = result.vry;
                    let ztd = result.vrz;

                    // Calculate target states at desired intercept time
                    let x0 = StateVector::new(
                        xt / 3280.0, yt / 3280.0, zt / 3280.0,
                        xtd / 3280.0, ytd / 3280.0, ztd / 3280.0,
                    );
                    let x1 = kepler1(&x0, 0.0, tf);
                    let xtf = x1.x * 3280.0;
                    let ytf = x1.y * 3280.0;
                    let ztf = x1.z * 3280.0;

                    let altfkm = ((xt * xt + ytf * ytf + ztf * ztf).sqrt() - a) / 3280.0;
                    if altfkm >= 50.0 {
                        // Calculate missile velocity required
                        let tgolamm = tf - tlaunch;
                        let result_m = lambert3d(xm, ym, zm, tgolamm, xtf, ytf, ztf, switchm);
                        let xmd = result_m.vrx;
                        let ymd = result_m.vry;
                        let zmd = result_m.vrz;

                        let vbom = (xmd * xmd + ymd * ymd + zmd * zmd).sqrt() / 3280.0;
                        let vbot = (xtd * xtd + ytd * ytd + ztd * ztd).sqrt() / 3280.0;

                        if vbom < vbolim && vbot < 7.5 {
                            array_tf.push(tf);
                            array_tlaunch.push(tlaunch);
                            array_xlongt.push(xlongtdeg * 111.0);
                            array_xlatt.push(xlattdeg * 111.0);
                            array_vbom.push(vbom);
                            array_vbot.push(vbot);
                        }
                    }

                    tf += 50.0;
                }
                tlaunch += 100.0;
            }
            xlattdeg += 2.5;
        }
        xlongtdeg += 5.0;
    }

    Results {
        tf: array_tf,
        tlaunch: array_tlaunch,
        xlongt: array_xlongt,
        xlatt: array_xlatt,
        vbom: array_vbom,
        vbot: array_vbot,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c44l4_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.tf.clone(),
        results.tlaunch.clone(),
        results.xlongt.clone(),
        results.xlatt.clone(),
        results.vbom.clone(),
        results.vbot.clone(),
    ])?;

    println!("C44L4: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c44l4_runs() {
        let results = run();
        assert!(true);
    }
}
