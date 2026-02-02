//! Chapter 17, Lesson 2: Lambert Solver Direct Call
//!
//! Direct call to lambertpz for orbit transfer.
//!
//! NOTE: The MATLAB code does NOT produce datfil output - it only prints VRXM and
//! VRYM values. This is a demonstration of the Lambert solver, not a simulation.
//! MATLAB BUG: Calls `lambertpz` but the function file is `lamberpz.m` and the
//! function is defined as `lambert` not `lambertpz`.

use crate::utils::lambert2d::lambert2d;

pub struct Results {
    pub vrxm: f64,
    pub vrym: f64,
}

/// Run the C17L2 simulation
pub fn run() -> Results {
    let xlongmdeg: f64 = 45.0;
    let xlongtdeg: f64 = 90.0;
    let altnmt: f64 = 0.0;
    let altnmm: f64 = 0.0;
    let tf: f64 = 1000.0;
    let pi: f64 = std::f64::consts::PI;
    let degrad = 360.0 / (2.0 * pi);
    let a: f64 = 2.0926e7;

    let altt = altnmt * 6076.0;
    let altm = altnmm * 6076.0;
    let xlongm = xlongmdeg / degrad;
    let xlongt = xlongtdeg / degrad;
    let xm = (a + altm) * xlongm.cos();
    let ym = (a + altm) * xlongm.sin();
    let xt = (a + altt) * xlongt.cos();
    let yt = (a + altt) * xlongt.sin();

    let result = lambert2d(xm, ym, tf, xt, yt, xlongm, xlongt);

    println!("VRXM = {:.6e}", result.vrx);
    println!("VRYM = {:.6e}", result.vry);

    Results {
        vrxm: result.vrx,
        vrym: result.vry,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let _data_file = format!("{}/c17l2_datfil.txt", output_dir);
    // No array output for this simulation

    println!("C17L2: Simulation finished");
    println!("  VRXM = {:.6e}", results.vrxm);
    println!("  VRYM = {:.6e}", results.vrym);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c17l2_runs() {
        let _results = run();
    }
}
