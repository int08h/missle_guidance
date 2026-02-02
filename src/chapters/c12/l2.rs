//! Chapter 12, Lesson 2: Extended Kalman Filter for Beta Estimation
//!
//! EKF estimating ballistic coefficient (beta) of reentry vehicle.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;
use rand_distr::{Distribution, Normal};

pub struct Results {
    pub time: Vec<f64>,
    pub x: Vec<f64>,
    pub xh: Vec<f64>,
    pub xd: Vec<f64>,
    pub xdh: Vec<f64>,
    pub beta: Vec<f64>,
    pub betah: Vec<f64>,
    pub errx: Vec<f64>,
    pub sp11: Vec<f64>,
    pub errxd: Vec<f64>,
    pub sp22: Vec<f64>,
    pub errbetainv: Vec<f64>,
    pub sp33: Vec<f64>,
}

/// Project state forward using Euler integration
fn project_c12l2(ts: f64, xp: f64, xdp: f64, betainvp: f64, hp: f64) -> (f64, f64) {
    let mut t: f64 = 0.0;
    let mut x = xp;
    let mut xd = xdp;
    let betainv = betainvp;
    let h = hp;

    while t <= (ts - 0.0001) {
        let xdd = 0.0034 * 32.2 * xd * xd * betainv * (-x / 22000.0).exp() / 2.0 - 32.2;
        xd += h * xdd;
        x += h * xd;
        t += h;
    }

    (x, xd)
}

/// Run the C12L2 simulation
pub fn run() -> Results {
    let g: f64 = 32.2;
    let signoise: f64 = 25.0;
    let mut x: f64 = 200000.0;
    let mut xd: f64 = -6000.0;
    let beta: f64 = 500.0;
    let mut xh: f64 = 200025.0;
    let mut xdh: f64 = -6150.0;
    let mut betainvh: f64 = 1.0 / 800.0;
    let betainv: f64 = 1.0 / beta;
    let _order: usize = 3;
    let ts: f64 = 0.1;
    let tf: f64 = 30.0;
    let _phis: f64 = 0.0;
    let h: f64 = 0.001;
    let hp: f64 = 0.001;

    let normal = Normal::new(0.0, 1.0).unwrap();
    let mut rng = rand::thread_rng();

    // Covariance matrix P (stored as elements)
    let mut p11 = signoise * signoise;
    let mut p12: f64 = 0.0;
    let mut p13: f64 = 0.0;
    let mut p21: f64 = 0.0;
    let mut p22: f64 = 20000.0;
    let mut p23: f64 = 0.0;
    let mut p31: f64 = 0.0;
    let mut p32: f64 = 0.0;
    let mut p33: f64 = (betainv - betainvh).powi(2);

    let rmat = signoise * signoise;

    let mut t: f64 = 0.0;
    let mut s: f64 = 0.0;

    let mut array_t = Vec::new();
    let mut array_x = Vec::new();
    let mut array_xh = Vec::new();
    let mut array_xd = Vec::new();
    let mut array_xdh = Vec::new();
    let mut array_beta = Vec::new();
    let mut array_betah = Vec::new();
    let mut array_errx = Vec::new();
    let mut array_sp11 = Vec::new();
    let mut array_errxd = Vec::new();
    let mut array_sp22 = Vec::new();
    let mut array_errbetainv = Vec::new();
    let mut array_sp33 = Vec::new();

    while t <= tf {
        let x_old = x;
        let xd_old = xd;

        // First derivative
        let xdd = 0.0034 * g * xd * xd * (-x / 22000.0).exp() / (2.0 * beta) - g;
        x += h * xd;
        xd += h * xdd;
        t += h;

        // Second derivative
        let xdd = 0.0034 * g * xd * xd * (-x / 22000.0).exp() / (2.0 * beta) - g;
        x = 0.5 * (x_old + x + h * xd);
        xd = 0.5 * (xd_old + xd + h * xdd);

        s += h;
        if s >= (ts - 0.00001) {
            s = 0.0;

            // Compute F matrix elements
            let rhoh = 0.0034 * (-xh / 22000.0).exp();
            let f21 = -g * rhoh * xdh * xdh * betainvh / 44000.0;
            let f22 = rhoh * g * xdh * betainvh;
            let f23 = 0.5 * rhoh * xdh * xdh * g;

            // PHI = I + F*TS (first-order approximation)
            let phi11 = 1.0;
            let phi12 = ts;
            let phi13 = 0.0;
            let phi21 = f21 * ts;
            let phi22 = 1.0 + f22 * ts;
            let phi23 = f23 * ts;
            let phi31 = 0.0;
            let phi32 = 0.0;
            let phi33 = 1.0;

            // Q matrix (process noise covariance)
            let phis: f64 = 0.0;
            let q11 = 0.0;
            let q12 = 0.0;
            let q13 = 0.0;
            let q21 = 0.0;
            let q22 = f23 * f23 * phis * ts * ts * ts / 3.0;
            let q23 = f23 * phis * ts * ts / 2.0;
            let q31 = 0.0;
            let q32 = f23 * phis * ts * ts / 2.0;
            let q33 = phis * ts;

            // M = PHI * P * PHI' + Q
            // First compute PHI * P
            let pp11 = phi11 * p11 + phi12 * p21 + phi13 * p31;
            let pp12 = phi11 * p12 + phi12 * p22 + phi13 * p32;
            let pp13 = phi11 * p13 + phi12 * p23 + phi13 * p33;
            let pp21 = phi21 * p11 + phi22 * p21 + phi23 * p31;
            let pp22 = phi21 * p12 + phi22 * p22 + phi23 * p32;
            let pp23 = phi21 * p13 + phi22 * p23 + phi23 * p33;
            let pp31 = phi31 * p11 + phi32 * p21 + phi33 * p31;
            let pp32 = phi31 * p12 + phi32 * p22 + phi33 * p32;
            let pp33 = phi31 * p13 + phi32 * p23 + phi33 * p33;

            // Then (PHI * P) * PHI'
            let m11 = pp11 * phi11 + pp12 * phi12 + pp13 * phi13 + q11;
            let m12 = pp11 * phi21 + pp12 * phi22 + pp13 * phi23 + q12;
            let m13 = pp11 * phi31 + pp12 * phi32 + pp13 * phi33 + q13;
            let m21 = pp21 * phi11 + pp22 * phi12 + pp23 * phi13 + q21;
            let m22 = pp21 * phi21 + pp22 * phi22 + pp23 * phi23 + q22;
            let m23 = pp21 * phi31 + pp22 * phi32 + pp23 * phi33 + q23;
            let m31 = pp31 * phi11 + pp32 * phi12 + pp33 * phi13 + q31;
            let m32 = pp31 * phi21 + pp32 * phi22 + pp33 * phi23 + q32;
            let m33 = pp31 * phi31 + pp32 * phi32 + pp33 * phi33 + q33;

            // K = M * H' / (H * M * H' + R)
            // H = [1 0 0], so H*M*H' = M(1,1)
            let k1 = m11 / (m11 + rmat);
            let k2 = m21 / (m11 + rmat);
            let k3 = m31 / (m11 + rmat);

            // P = (I - K*H) * M
            p11 = (1.0 - k1) * m11;
            p12 = (1.0 - k1) * m12;
            p13 = (1.0 - k1) * m13;
            p21 = -k2 * m11 + m21;
            p22 = -k2 * m12 + m22;
            p23 = -k2 * m13 + m23;
            p31 = -k3 * m11 + m31;
            p32 = -k3 * m12 + m32;
            p33 = -k3 * m13 + m33;

            let xnoise = signoise * normal.sample(&mut rng);

            // Project estimate forward
            let (xb, xdb) = project_c12l2(ts, xh, xdh, betainvh, hp);
            let res = x + xnoise - xb;

            xh = xb + k1 * res;
            xdh = xdb + k2 * res;
            betainvh += k3 * res;

            let errx = x - xh;
            let sp11_val = p11.sqrt();
            let errxd = xd - xdh;
            let sp22_val = p22.sqrt();
            let errbetainv_val = betainv - betainvh;
            let sp33_val = p33.sqrt();
            let betah_val = 1.0 / betainvh;

            array_t.push(t);
            array_x.push(x);
            array_xh.push(xh);
            array_xd.push(xd);
            array_xdh.push(xdh);
            array_beta.push(beta);
            array_betah.push(betah_val);
            array_errx.push(errx);
            array_sp11.push(sp11_val);
            array_errxd.push(errxd);
            array_sp22.push(sp22_val);
            array_errbetainv.push(errbetainv_val);
            array_sp33.push(sp33_val);
        }
    }

    Results {
        time: array_t,
        x: array_x,
        xh: array_xh,
        xd: array_xd,
        xdh: array_xdh,
        beta: array_beta,
        betah: array_betah,
        errx: array_errx,
        sp11: array_sp11,
        errxd: array_errxd,
        sp22: array_sp22,
        errbetainv: array_errbetainv,
        sp33: array_sp33,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c12l2_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.x.clone(),
        results.xh.clone(),
        results.xd.clone(),
        results.xdh.clone(),
        results.beta.clone(),
        results.betah.clone(),
    ])?;

    let cov_file = format!("{}/c12l2_covfil.txt", output_dir);
    let sp11p: Vec<f64> = results.sp11.iter().map(|v| -v).collect();
    let sp22p: Vec<f64> = results.sp22.iter().map(|v| -v).collect();
    let sp33p: Vec<f64> = results.sp33.iter().map(|v| -v).collect();
    save_data(&cov_file, &[
        results.time.clone(),
        results.errx.clone(),
        results.sp11.clone(),
        sp11p,
        results.errxd.clone(),
        results.sp22.clone(),
        sp22p,
        results.errbetainv.clone(),
        results.sp33.clone(),
        sp33p,
    ])?;

    // Altitude error plot
    let plot_file = format!("{}/c12l2_errx.png", output_dir);
    let config = PlotConfig::new("Error in Estimate of Altitude")
        .with_labels("Time (Sec)", "Error in Estimate of Altitude (Ft)")
        .with_y_range(-25.0, 25.0);
    let series = vec![
        Series::new(results.time.clone(), results.errx.clone())
            .with_label("Error")
            .with_color(plotters::prelude::BLUE),
        Series::new(results.time.clone(), results.sp11.clone())
            .with_label("+Sigma")
            .with_color(plotters::prelude::RED),
    ];
    line_plot(&plot_file, &config, &series).ok();

    println!("C12L2: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c12l2_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
