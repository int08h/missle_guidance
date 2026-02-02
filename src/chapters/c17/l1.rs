//! Chapter 17, Lesson 1: Lambert Orbit Solver
//!
//! Use olambert to find initial velocity for orbit transfer.
//! Outputs iteration history (ArrayVRX, ArrayVRY, ArrayTF).

use crate::save_data;
use crate::utils::olambert::olambert;

use std::f64::consts::PI;

pub struct Results {
    pub vrx: Vec<f64>,
    pub vry: Vec<f64>,
    pub tf: Vec<f64>,
}

/// Run the C17L1 simulation
pub fn run() -> Results {
    let xlongmdeg: f64 = 45.0;
    let xlongtdeg: f64 = 90.0;
    let altnmt: f64 = 0.0;
    let altnmm: f64 = 0.0;
    let tf_val: f64 = 1000.0;
    let degrad = 360.0 / (2.0 * PI);
    let a: f64 = 2.0926e7;
    let gm: f64 = 1.4077e16;

    let altt = altnmt * 6076.0;
    let altm = altnmm * 6076.0;
    let xlongm = xlongmdeg / degrad;
    let xlongt = xlongtdeg / degrad;
    let xm = (a + altm) * xlongm.cos();
    let ym = (a + altm) * xlongm.sin();
    let xt = (a + altt) * xlongt.cos();
    let yt = (a + altt) * xlongt.sin();

    // Use olambert solver (matches MATLAB olambert.m)
    let result = olambert(xm, ym, tf_val, xt, yt, xlongm, xlongt, a, gm);

    println!("The final iteration");
    println!("count = {}", result.count);
    println!("VRXM = {:.6e}", result.vrx);
    println!("VRYM = {:.6e}", result.vry);

    // Return iteration history arrays (matching MATLAB output)
    Results {
        vrx: result.array_vrx,
        vry: result.array_vry,
        tf: result.array_tf,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c17l1_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.vrx.clone(),
        results.vry.clone(),
        results.tf.clone(),
    ])?;

    println!("C17L1: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c17l1_runs() {
        let results = run();
        assert!(!results.vrx.is_empty());
    }
}
