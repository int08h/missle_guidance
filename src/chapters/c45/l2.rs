//! Chapter 45, Lesson 2: TMD IRBM Trajectory
//!
//! Theater missile defense: IRBM trajectory simulation with Lambert guidance.

use crate::save_data;
use crate::utils::{lambert3d, distance3dkm, EARTH_RADIUS_FT, GM_FT};

pub struct Results {
    pub time: Vec<f64>,
    pub dist_km: Vec<f64>,
    pub alt_km: Vec<f64>,
}

/// Run the C45L2 simulation
pub fn run() -> Results {
    let rdeskm: f64 = 2000.0;
    let itgt: i32 = 1;  // 1=IRBM, 2=ICBM
    let tloft: f64 = 200.0;
    let tupt: f64 = 15.0;
    let switch: i32 = 0;

    let tpz = if itgt == 1 { 180.0 } else { 240.0 };

    let a = EARTH_RADIUS_FT;
    let gm = GM_FT;

    let gamdeg: f64 = 89.99;
    let h: f64 = 0.01;
    let mut t: f64 = 0.0;
    let mut s: f64 = 0.0;
    let alt: f64 = 0.0;
    let angdeg: f64 = 0.0;
    let ang = angdeg / 57.3;

    let mut x = (a + alt) * ang.cos();
    let mut y = (a + alt) * ang.sin();
    let z: f64 = 0.0;
    let mut alt = (x * x + y * y).sqrt() - a;

    let xfirst = x;
    let yfirst = y;
    let zfirst = z;

    let mut x1 = (1.5708 - gamdeg / 57.3 + ang).cos();
    let mut y1 = (1.5708 - gamdeg / 57.3 + ang).sin();

    let mut axt: f64 = 0.0;
    let mut ayt: f64 = 0.0;

    let xlongtdeg = 57.3 * rdeskm * 3280.0 / a;
    let mut tf = 252.0 + 0.223 * rdeskm - 5.44e-6 * rdeskm * rdeskm;
    tf += tloft;

    let xlongt = xlongtdeg / 57.3;
    let xf = a * xlongt.cos();
    let yf = a * xlongt.sin();
    let zf: f64 = 0.0;

    let mut qboost: bool = true;

    let mut array_t = Vec::new();
    let mut array_distkm = Vec::new();
    let mut array_altkm = Vec::new();

    while alt > -1.0 {
        let xold = x;
        let yold = y;
        let x1old = x1;
        let y1old = y1;

        // Thrust model
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

        let at = 32.2 * trst / wgt;
        let tembot = (x * x + y * y).powf(1.5);
        let x1d = -gm * x / tembot + axt;
        let y1d = -gm * y / tembot + ayt;
        alt = (x * x + y * y).sqrt() - a;

        // Euler step
        x += h * x1;
        y += h * y1;
        x1 += h * x1d;
        y1 += h * y1d;
        t += h;

        // RK2 averaging
        x = (xold + x) / 2.0 + 0.5 * h * x1;
        y = (yold + y) / 2.0 + 0.5 * h * y1;
        x1 = (x1old + x1) / 2.0 + 0.5 * h * x1d;
        y1 = (y1old + y1) / 2.0 + 0.5 * h * y1d;

        s += h;

        // Lambert guidance
        if qboost {
            let tgolam = tf - t;
            let result = lambert3d(x, y, z, tgolam, xf, yf, zf, switch);
            let vrx = result.vrx;
            let vry = result.vry;

            let delx = vrx - x1;
            let dely = vry - y1;
            let del = (delx * delx + dely * dely).sqrt();

            if t < tpz && del > 500.0 {
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

            if t < tupt {
                let rtmag = (x * x + y * y).sqrt();
                axt = at * x / rtmag;
                ayt = at * y / rtmag;
            }
        }

        if s >= 0.99999 {
            s = 0.0;
            let distkm = distance3dkm(x, y, z, xfirst, yfirst, zfirst);
            let altkm = ((x * x + y * y).sqrt() - a) / 3280.0;
            let _velk = (x1 * x1 + y1 * y1).sqrt() / 3280.0;

            array_t.push(t);
            array_distkm.push(distkm);
            array_altkm.push(altkm);
        }
    }

    Results {
        time: array_t,
        dist_km: array_distkm,
        alt_km: array_altkm,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c45l2_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.dist_km.clone(),
        results.alt_km.clone(),
    ])?;

    println!("C45L2: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c45l2_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
