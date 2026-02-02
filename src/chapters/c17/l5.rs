//! Chapter 17, Lesson 5: Ballistic Range Calculation
//!
//! Iterative calculation of final range for given initial conditions.

use crate::save_data;

pub struct Results {
    pub icount: Vec<f64>,
    pub phideg: Vec<f64>,
    pub rf: Vec<f64>,
    pub tf: Vec<f64>,
}

/// Run the C17L5 simulation
pub fn run() -> Results {
    let gamdeg: f64 = 30.0;
    let v: f64 = 24000.0;
    let tfdes: f64 = 2200.0;
    let pi: f64 = std::f64::consts::PI;
    let degrad = 360.0 / (2.0 * pi);
    let xlongmdeg: f64 = 0.0;
    let a: f64 = 2.0926e7;
    let gm: f64 = 1.4077e16;
    let gam = gamdeg / degrad;
    let _xlongm = xlongmdeg / degrad;
    let altm: f64 = 0.0;
    let xm = (a + altm) * (xlongmdeg / degrad).cos();
    let ym = (a + altm) * (xlongmdeg / degrad).sin();
    let r0 = (xm * xm + ym * ym).sqrt();

    let mut phimax = pi / 2.0;
    let mut phimin: f64 = 0.0;
    let mut phi = 45.0 / degrad;
    let mut tf: f64 = 100000.0;

    let mut array_icount = Vec::new();
    let mut array_phideg = Vec::new();
    let mut array_rf = Vec::new();
    let mut array_tf = Vec::new();

    let mut icount: usize = 0;
    let mut phiold: f64 = 0.0;
    let mut told: f64 = 0.0;

    while (tfdes - tf).abs() > 0.00000001 * tfdes {
        let phideg_val = phi * degrad;
        let top = v * v * r0 * r0 * gam.cos() * gam.cos();
        let bot = gm * (1.0 - phi.cos()) + r0 * v * v * gam.cos() * (phi + gam).cos();
        let rf_val = top / bot;
        let xlam = r0 * v * v / gm;

        let top1 = gam.tan() * (1.0 - phi.cos()) + (1.0 - xlam) * phi.sin();
        let bot1p = (1.0 - phi.cos()) / (xlam * gam.cos() * gam.cos());
        let bot1 = (2.0 - xlam) * (bot1p + (gam + phi).cos() / gam.cos());
        let top2 = 2.0 * gam.cos();
        let bot2 = xlam * (2.0 / xlam - 1.0).powf(1.5);
        let top3 = (2.0 / xlam - 1.0).sqrt();
        let bot3 = gam.cos() / (phi / 2.0).tan() - gam.sin();
        let temp = (top2 / bot2) * top3.atan2(bot3);
        tf = r0 * (top1 / bot1 + temp) / (v * gam.cos());

        icount += 1;

        if tf > tfdes {
            phimax = phi;
        } else {
            phimin = phi;
        }

        let xnext = if icount == 1 {
            (phimax + phimin) / 2.0
        } else {
            let mut xnext = phi + (phi - phiold) * (tfdes - tf) / (tf - told);
            if xnext > phimax || xnext < phimin {
                xnext = (phimax + phimin) / 2.0;
            }
            xnext
        };

        phiold = phi;
        told = tf;
        phi = xnext;

        array_icount.push(icount as f64);
        array_phideg.push(phideg_val);
        array_rf.push(rf_val);
        array_tf.push(tf);
    }

    let xf = array_rf.last().unwrap_or(&0.0) * phi.cos();
    let yf = array_rf.last().unwrap_or(&0.0) * phi.sin();
    println!("XF = {:.6e}", xf);
    println!("YF = {:.6e}", yf);

    Results {
        icount: array_icount,
        phideg: array_phideg,
        rf: array_rf,
        tf: array_tf,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c17l5_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.icount.clone(),
        results.phideg.clone(),
        results.rf.clone(),
        results.tf.clone(),
    ])?;

    println!("C17L5: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c17l5_runs() {
        let results = run();
        assert!(!results.icount.is_empty());
    }
}
