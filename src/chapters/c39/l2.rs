//! Chapter 39, Lesson 2: Optimal Guidance Miss Analysis
//!
//! Miss distance analysis for optimal guidance with autopilot dynamics.

use crate::save_data;

pub struct Results {
    pub tf: Vec<f64>,
    pub miss: Vec<f64>,
}

/// Run the C39L2 simulation
pub fn run() -> Results {
    let vc: f64 = 4000.0;
    let gam: f64 = 0.0001;
    let apn: i32 = 2;
    let xnt: f64 = 161.0;
    let tau: f64 = 0.5;
    let wz: f64 = 10.0;
    let xnp: f64 = 3.0;
    let xnclim: f64 = 9999999.0;
    let h: f64 = 0.01;

    let mut array_tf = Vec::new();
    let mut array_y = Vec::new();

    let mut tf = 0.1;
    while tf <= 10.0 {
        let mut y: f64 = 0.0;
        let mut yd: f64 = 0.0;
        let mut e: f64 = 0.0;
        let mut t: f64 = 0.0;
        let mut xnl: f64 = 0.0;

        while t <= tf - 1e-5 {
            let yold = y;
            let ydold = yd;
            let eold = e;

            // First derivative evaluation
            let tgo = tf - t + 0.00001;
            let xlamd = (y + yd * tgo) / (vc * tgo * tgo);

            let mut xnc = if apn == 0 {
                xnp * (y + yd * tgo) / (tgo * tgo)
            } else if apn == 1 {
                let x = tgo / tau;
                let top = 6.0 * x * x * ((-x).exp() - 1.0 + x);
                let bot1 = 2.0 * x * x * x + 3.0 + 6.0 * x - 6.0 * x * x;
                let bot2 = -12.0 * x * (-x).exp() - 3.0 * (-2.0 * x).exp();
                let xnpp = top / (0.0001 + bot1 + bot2);
                let xnew = xnpp * xnl * ((-x).exp() + x - 1.0) / (x * x);
                xnpp * vc * xlamd + 0.5 * xnpp * xnt - xnew
            } else {
                let xs = tgo / tau;
                let temp1 = tgo * tgo * tau * ((-xs).exp() - 1.0 + xs);
                let top = -(tgo.powi(3)) / (tau * wz) + (1.0 + 1.0 / (tau * wz)) * temp1;
                let temp2 = 0.5 * (1.0 - 3.0 / (tau * wz)) + xs * (1.0 + 1.0 / (tau * wz)) - xs * xs;
                let temp3 = -2.0 * xs * (-xs).exp();
                let temp4 = 2.0 * (-xs).exp() / (tau * wz) - 0.5 * (-2.0 * xs).exp() * (1.0 + 1.0 / (tau * wz));
                let bot = gam + tgo.powi(3) / 3.0 + (1.0 + 1.0 / (tau * wz)) * tau.powi(3) * (temp2 + temp3 + temp4);
                let xnpp = top / bot;
                let c1th = xnpp / (tgo * tgo);
                let c2th = xnpp / tgo;
                let c3th = 0.5 * xnpp;
                let c4th = -xnpp * ((-xs).exp() + xs - 1.0) / (xs * xs);
                c1th * y + c2th * yd + c3th * xnt + c4th * e
            };

            xnc = xnc.clamp(-xnclim, xnclim);

            let ed = (xnc - e) / tau;
            xnl = e - ed / wz;
            let ydd = xnt - xnl;

            // Euler step
            y += h * yd;
            yd += h * ydd;
            e += h * ed;
            t += h;

            // Second evaluation for RK2
            let tgo = tf - t + 0.00001;
            let xlamd = (y + yd * tgo) / (vc * tgo * tgo);

            let mut xnc = if apn == 0 {
                xnp * (y + yd * tgo) / (tgo * tgo)
            } else if apn == 1 {
                let x = tgo / tau;
                let top = 6.0 * x * x * ((-x).exp() - 1.0 + x);
                let bot1 = 2.0 * x * x * x + 3.0 + 6.0 * x - 6.0 * x * x;
                let bot2 = -12.0 * x * (-x).exp() - 3.0 * (-2.0 * x).exp();
                let xnpp = top / (0.0001 + bot1 + bot2);
                let xnew = xnpp * xnl * ((-x).exp() + x - 1.0) / (x * x);
                xnpp * vc * xlamd + 0.5 * xnpp * xnt - xnew
            } else {
                let xs = tgo / tau;
                let temp1 = tgo * tgo * tau * ((-xs).exp() - 1.0 + xs);
                let top = -(tgo.powi(3)) / (tau * wz) + (1.0 + 1.0 / (tau * wz)) * temp1;
                let temp2 = 0.5 * (1.0 - 3.0 / (tau * wz)) + xs * (1.0 + 1.0 / (tau * wz)) - xs * xs;
                let temp3 = -2.0 * xs * (-xs).exp();
                let temp4 = 2.0 * (-xs).exp() / (tau * wz) - 0.5 * (-2.0 * xs).exp() * (1.0 + 1.0 / (tau * wz));
                let bot = gam + tgo.powi(3) / 3.0 + (1.0 + 1.0 / (tau * wz)) * tau.powi(3) * (temp2 + temp3 + temp4);
                let xnpp = top / bot;
                let c1th = xnpp / (tgo * tgo);
                let c2th = xnpp / tgo;
                let c3th = 0.5 * xnpp;
                let c4th = -xnpp * ((-xs).exp() + xs - 1.0) / (xs * xs);
                c1th * y + c2th * yd + c3th * xnt + c4th * e
            };

            xnc = xnc.clamp(-xnclim, xnclim);

            let ed = (xnc - e) / tau;
            xnl = e - ed / wz;
            let ydd = xnt - xnl;

            // RK2 averaging
            y = 0.5 * (yold + y + h * yd);
            yd = 0.5 * (ydold + yd + h * ydd);
            e = 0.5 * (eold + e + h * ed);
        }

        array_tf.push(tf);
        array_y.push(y);
        tf += 0.1;
    }

    Results {
        tf: array_tf,
        miss: array_y,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data_file = format!("{}/c39l2_datfil.txt", output_dir);
    save_data(&data_file, &[
        results.tf.clone(),
        results.miss.clone(),
    ])?;

    println!("C39L2: Simulation finished");
    println!("  Data saved to: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c39l2_runs() {
        let results = run();
        assert!(!results.tf.is_empty());
    }
}
