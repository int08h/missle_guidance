//! Chapter 20, Listing 9: Adjoint with Second-Order Lag (Sweep THOM with QSWITCH)
//!
//! Simulates adjoint response sweeping through different homing times (THOM).
//! Uses a second-order lag filter with QSWITCH logic for displacement.

use crate::save_data;

pub struct Results {
    pub thom: Vec<f64>,
    pub y: Vec<f64>,
}

/// Run the C20L9 simulation
pub fn run() -> Results {
    let vc: f64 = 4000.0;
    let xnt: f64 = 0.0;
    let displace: f64 = 200.0;
    let _vm: f64 = 3000.0;
    let tau: f64 = 1.0;
    let xnp: f64 = 3.0;
    let xnclim: f64 = 99999999.0;
    let tf: f64 = 10.0;

    let mut array_thom = Vec::new();
    let mut array_y = Vec::new();

    // Sweep THOM from 0.1 to 10.0 in steps of 0.1
    let mut thom_val: f64 = 0.1;
    while thom_val <= 10.0 + 0.00001 {
        let mut qswitch = false;
        let mut y: f64 = 0.0;
        let mut yd: f64 = 0.0;
        let mut xnl: f64 = 0.0;
        let mut d: f64 = 0.0;
        let mut elamdh: f64 = 0.0;
        let mut x4: f64 = 0.0;
        let mut x5: f64 = 0.0;
        let mut t: f64 = 0.0;
        let h: f64 = 0.01;

        while t <= tf - 0.0001 {
            let mut tgo = tf - t + 0.00001;

            // Check for QSWITCH (displacement at specified TGO)
            if tgo <= thom_val && !qswitch {
                qswitch = true;
                y += displace;
                let xlam = y / (vc * tgo);
                d = xlam;
            }

            let yold = y;
            let ydold = yd;
            let xnlold = xnl;
            let dold = d;
            let elamdhold = elamdh;
            let x4old = x4;
            let x5old = x5;

            // First derivative evaluation
            tgo = tf - t + 0.00001;
            let mut xlam = y / (vc * tgo);
            let mut dd = 5.0 * (xlam - d) / tau;
            let mut elamdhd = 5.0 * (dd - elamdh) / tau;
            let mut xnc = xnp * vc * elamdh;
            if xnc > xnclim {
                xnc = xnclim;
            }
            if xnc < -xnclim {
                xnc = -xnclim;
            }
            let mut x4d = 5.0 * (xnc - x4) / tau;
            let mut x5d = 5.0 * (x4 - x5) / tau;
            let mut xnld = 5.0 * (x5 - xnl) / tau;
            let mut ydd = xnt - xnl;

            // Euler step
            y += h * yd;
            yd += h * ydd;
            xnl += h * xnld;
            elamdh += h * elamdhd;
            d += h * dd;
            x4 += h * x4d;
            x5 += h * x5d;
            t += h;

            // Second derivative evaluation
            tgo = tf - t + 0.00001;
            xlam = y / (vc * tgo);
            dd = 5.0 * (xlam - d) / tau;
            elamdhd = 5.0 * (dd - elamdh) / tau;
            xnc = xnp * vc * elamdh;
            if xnc > xnclim {
                xnc = xnclim;
            }
            if xnc < -xnclim {
                xnc = -xnclim;
            }
            x4d = 5.0 * (xnc - x4) / tau;
            x5d = 5.0 * (x4 - x5) / tau;
            xnld = 5.0 * (x5 - xnl) / tau;
            ydd = xnt - xnl;

            // RK2 averaging
            y = 0.5 * (yold + y + h * yd);
            yd = 0.5 * (ydold + yd + h * ydd);
            xnl = 0.5 * (xnlold + xnl + h * xnld);
            d = 0.5 * (dold + d + h * dd);
            elamdh = 0.5 * (elamdhold + elamdh + h * elamdhd);
            x4 = 0.5 * (x4old + x4 + h * x4d);
            x5 = 0.5 * (x5old + x5 + h * x5d);
        }

        array_thom.push(thom_val);
        array_y.push(y);

        thom_val += 0.1;
    }

    Results {
        thom: array_thom,
        y: array_y,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data: Vec<Vec<f64>> = vec![
        results.thom.clone(),
        results.y.clone(),
    ];
    let data_file = format!("{}/c20l9_datfil.txt", output_dir);
    save_data(&data_file, &data)?;

    println!("C20L9: Simulation complete");
    println!("  Output: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c20l9_runs() {
        let results = run();
        assert!(!results.thom.is_empty());
    }
}
