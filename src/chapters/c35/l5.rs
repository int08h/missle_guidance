//! Chapter 35, Lesson 5: Adjoint Model for Autopilot
//!
//! Computes optimal guidance gains using adjoint method
//! with third-order autopilot dynamics.

use crate::plotting::{PlotConfig, Series, line_plot};
use crate::save_data;

pub struct Results {
    pub time: Vec<f64>,
    pub c4: Vec<f64>,
    pub c5: Vec<f64>,
    pub c6: Vec<f64>,
    pub xnp: Vec<f64>,
}

/// Run the C35L5 simulation
pub fn run() -> Results {
    let tau: f64 = 1.0;
    let gam: f64 = 0.00001;
    let wz: f64 = 5.0;
    let w: f64 = 20.0;
    let z: f64 = 0.7;

    let mut t: f64 = 0.0;
    let h: f64 = 0.01;
    let mut s: f64 = 0.0;

    // State variables
    let mut x: f64 = 0.0;
    let mut e: f64 = 0.0;
    let mut ed: f64 = 0.0;
    let mut edd: f64 = 0.0;
    let mut e1: f64 = 0.0;
    let mut e2: f64 = 0.0;
    let mut e2d: f64 = 0.0;
    let mut e3: f64 = 0.0;
    let mut e3d: f64 = 0.0;
    let mut e4: f64 = 0.0;
    let mut e4d: f64;
    let mut e5: f64 = 0.0;
    let mut e5d: f64 = 0.0;
    let mut e5dd: f64 = 0.0;

    let mut array_t = Vec::new();
    let mut array_c4 = Vec::new();
    let mut array_c5 = Vec::new();
    let mut array_c6 = Vec::new();
    let mut array_xnp = Vec::new();

    while t < 10.0 - 0.0001 {
        s += h;
        let xold = x;
        let eold = e;
        let edold = ed;
        let eddold = edd;
        let e1old = e1;
        let e2old = e2;
        let e2dold = e2d;
        let e3old = e3;
        let e3dold = e3d;
        let e4old = e4;
        let e5old = e5;
        let e5dold = e5d;
        let e5ddold = e5dd;

        // First derivative evaluation
        let eddd = w * w * (t - (1.0 / (w * w) + 2.0 * z * tau / w) * edd
            - (2.0 * z / w + tau) * ed - e) / tau;
        let xn = e - edd / (wz * wz);
        let xd = xn * xn;
        let mut d = x + gam;
        let mut pz = xn / d;
        let _xnp_val = t * t * pz;
        let e2dd = w * w * (-(1.0 / (w * w) + 2.0 * z * tau / w) * e2d
            - (2.0 * z / w + tau) * e2 - e1) / tau;
        let e1d = e2 + t;
        let _c4 = -(e1 - e2d / (wz * wz)) * pz;
        e4d = w * w * (-(1.0 / (w * w) + 2.0 * z * tau / w) * e4
            - (2.0 * z / w + tau) * e3d - e3) / tau;
        let e3dd = e4 + t;
        let _c5 = -(e3 - e4 / (wz * wz)) * pz;
        let e5ddd = t - w * w * ((1.0 / (w * w) + 2.0 * z * tau / w) * e5dd
            + (2.0 * z / w + tau) * e5d + e5) / tau;
        let _c6 = pz * (-e5 + e5dd / (wz * wz));

        // Euler step
        x += h * xd;
        e += h * ed;
        ed += h * edd;
        edd += h * eddd;
        e1 += h * e1d;
        e2 += h * e2d;
        e2d += h * e2dd;
        e3 += h * e3d;
        e3d += h * e3dd;
        e4 += h * e4d;
        e5 += h * e5d;
        e5d += h * e5dd;
        e5dd += h * e5ddd;
        t += h;

        // Second derivative for RK2
        // Note: In MATLAB, C4, C5, C6, XNP are recomputed here and these values are output
        let eddd = w * w * (t - (1.0 / (w * w) + 2.0 * z * tau / w) * edd
            - (2.0 * z / w + tau) * ed - e) / tau;
        let xn = e - edd / (wz * wz);
        let xd = xn * xn;
        d = x + gam;
        pz = xn / d;
        let xnp_val = t * t * pz;
        let e2dd = w * w * (-(1.0 / (w * w) + 2.0 * z * tau / w) * e2d
            - (2.0 * z / w + tau) * e2 - e1) / tau;
        let e1d = e2 + t;
        let c4 = -(e1 - e2d / (wz * wz)) * pz;
        e4d = w * w * (-(1.0 / (w * w) + 2.0 * z * tau / w) * e4
            - (2.0 * z / w + tau) * e3d - e3) / tau;
        let e3dd = e4 + t;
        let c5 = -(e3 - e4 / (wz * wz)) * pz;
        let e5ddd = t - w * w * ((1.0 / (w * w) + 2.0 * z * tau / w) * e5dd
            + (2.0 * z / w + tau) * e5d + e5) / tau;
        let c6 = pz * (-e5 + e5dd / (wz * wz));

        // RK2 averaging
        x = 0.5 * (xold + x + h * xd);
        e = 0.5 * (eold + e + h * ed);
        ed = 0.5 * (edold + ed + h * edd);
        edd = 0.5 * (eddold + edd + h * eddd);
        e1 = 0.5 * (e1old + e1 + h * e1d);
        e2 = 0.5 * (e2old + e2 + h * e2d);
        e2d = 0.5 * (e2dold + e2d + h * e2dd);
        e3 = 0.5 * (e3old + e3 + h * e3d);
        e3d = 0.5 * (e3dold + e3d + h * e3dd);
        e4 = 0.5 * (e4old + e4 + h * e4d);
        e5 = 0.5 * (e5old + e5 + h * e5d);
        e5d = 0.5 * (e5dold + e5d + h * e5dd);
        e5dd = 0.5 * (e5ddold + e5dd + h * e5ddd);

        if s >= 0.09999 {
            s = 0.0;
            array_t.push(t);
            array_c4.push(c4);
            array_c5.push(c5);
            array_c6.push(c6);
            array_xnp.push(xnp_val);
        }
    }

    Results {
        time: array_t,
        c4: array_c4,
        c5: array_c5,
        c6: array_c6,
        xnp: array_xnp,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c35l5_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.time.clone(),
        results.c4.clone(),
        results.c5.clone(),
        results.c6.clone(),
        results.xnp.clone(),
    ])?;

    let plot_file = format!("{}/c35l5_np.png", output_dir);
    let config = PlotConfig::new("Navigation Ratio NP")
        .with_labels("Time (s)", "NP");

    let series = vec![
        Series::new(results.time.clone(), results.xnp.clone())
            .with_color(plotters::prelude::BLUE),
    ];

    line_plot(&plot_file, &config, &series).ok();

    println!("C35L5: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c35l5_runs() {
        let results = run();
        assert!(!results.time.is_empty());
    }
}
