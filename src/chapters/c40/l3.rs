//! Chapter 40, Lesson 3: RK2 vs Kepler Comparison
//!
//! Compares RK2 numerical integration with Kepler analytical propagation.
//!
//! NOTE: The MATLAB code does NOT produce datfil output - it only prints comparison
//! values to the console. Calls KEPLER1 function for analytical propagation.

use crate::utils::{kepler1, StateVector, GM_FT};

#[allow(dead_code)]
pub struct Results {
    pub xtkm: f64,
    pub ytkm: f64,
    pub ztkm: f64,
    pub xtdkm: f64,
    pub ytdkm: f64,
    pub ztdkm: f64,
    pub x1_1: f64,
    pub x1_2: f64,
    pub x1_3: f64,
    pub x1_4: f64,
    pub x1_5: f64,
    pub x1_6: f64,
    pub errx: f64,
    pub erry: f64,
    pub errz: f64,
    pub errxd: f64,
    pub erryd: f64,
    pub errzd: f64,
}

/// Run the C40L3 simulation
pub fn run() -> Results {
    // Initial conditions from C28L2 output
    let mut xt: f64 = 14990432.9744621;
    let mut yt: f64 = 1952093.10305573;
    let mut zt: f64 = 14469752.1663352;
    let mut xtd: f64 = 996.773566434768;
    let mut ytd: f64 = -14954.8124604715;
    let mut ztd: f64 = 17528.1931768263;

    let tf: f64 = 1000.0;
    let gm = GM_FT;
    let h: f64 = 0.001;

    // Kepler propagation
    let t0: f64 = 0.0;
    let t1 = tf;
    let x0 = StateVector::new(
        xt / 3280.0,
        yt / 3280.0,
        zt / 3280.0,
        xtd / 3280.0,
        ytd / 3280.0,
        ztd / 3280.0,
    );
    let x1 = kepler1(&x0, t0, t1);

    // RK2 propagation
    let mut t = 0.0;
    let mut s = 0.0;

    while t <= tf {
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
        }
    }

    let xtkm = xt / 3280.0;
    let ytkm = yt / 3280.0;
    let ztkm = zt / 3280.0;
    let xtdkm = xtd / 3280.0;
    let ytdkm = ytd / 3280.0;
    let ztdkm = ztd / 3280.0;

    let errx = xtkm - x1.x;
    let erry = ytkm - x1.y;
    let errz = ztkm - x1.z;
    let errxd = xtdkm - x1.vx;
    let erryd = ytdkm - x1.vy;
    let errzd = ztdkm - x1.vz;

    Results {
        xtkm,
        ytkm,
        ztkm,
        xtdkm,
        ytdkm,
        ztdkm,
        x1_1: x1.x,
        x1_2: x1.y,
        x1_3: x1.z,
        x1_4: x1.vx,
        x1_5: x1.vy,
        x1_6: x1.vz,
        errx,
        erry,
        errz,
        errxd,
        erryd,
        errzd,
    }
}

pub fn run_and_save(_output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    // This listing doesn't produce datfil output in MATLAB
    // It just prints comparison values
    println!("C40L3: Simulation finished");
    println!("  XTKM (RK2):    {:.6}", results.xtkm);
    println!("  XTKM (Kepler): {:.6}", results.x1_1);
    println!("  YTKM (RK2):    {:.6}", results.ytkm);
    println!("  YTKM (Kepler): {:.6}", results.x1_2);
    println!("  ZTKM (RK2):    {:.6}", results.ztkm);
    println!("  ZTKM (Kepler): {:.6}", results.x1_3);
    println!("  ERRX: {:.6e}", results.errx);
    println!("  ERRY: {:.6e}", results.erry);
    println!("  ERRZ: {:.6e}", results.errz);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c40l3_runs() {
        let results = run();
        // Check that errors are small
        assert!(results.errx.abs() < 1.0);
    }
}
