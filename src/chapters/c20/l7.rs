//! Chapter 20, Listing 7: Adjoint with First-Order Lag (Heading Error)
//!
//! Simulates adjoint response including heading error contribution.
//! Computes total miss as sum of displacement and heading error terms.

use crate::save_data;

pub struct Results {
    pub tp: Vec<f64>,
    pub tot: Vec<f64>,
    pub xmhepz: Vec<f64>,
}

/// Run the C20L7 simulation
pub fn run() -> Results {
    let xnp: f64 = 3.0;
    let tau: f64 = 1.0;
    let tf: f64 = 5.0;
    let disp: f64 = 200.0;
    let vc: f64 = 4000.0;
    let vm: f64 = 3000.0;
    let he: f64 = -disp / vm;

    let mut _t: f64 = 0.0;
    let mut s: f64 = 0.0;
    let mut tp: f64 = _t + 0.00001;
    let mut x1: f64 = 0.0;
    let mut x2: f64 = 0.0;
    let mut x3: f64 = 1.0;
    let mut x4: f64 = 0.0;
    let h: f64 = 0.01;

    let mut array_tp = Vec::new();
    let mut array_tot = Vec::new();
    let mut array_xmhepz = Vec::new();

    while tp <= tf - 0.00001 {
        s += h;
        let x1old = x1;
        let x2old = x2;
        let x3old = x3;
        let x4old = x4;

        // First derivative evaluation
        let mut x1d = x2;
        let mut x2d = x3;
        let mut y1 = (x4 - xnp * vc * x2) / tau;
        let mut tgo = tp + 0.00001;
        let mut x3d = y1 / (vc * tgo);
        let mut x4d = -y1;

        // Euler step
        x1 += h * x1d;
        x2 += h * x2d;
        x3 += h * x3d;
        x4 += h * x4d;
        tp += h;

        // Second derivative evaluation
        x1d = x2;
        x2d = x3;
        y1 = (x4 - xnp * vc * x2) / tau;
        tgo = tp + 0.00001;
        x3d = y1 / (vc * tgo);
        x4d = -y1;

        // RK2 averaging
        x1 = (x1old + x1) / 2.0 + 0.5 * h * x1d;
        x2 = (x2old + x2) / 2.0 + 0.5 * h * x2d;
        x3 = (x3old + x3) / 2.0 + 0.5 * h * x3d;
        x4 = (x4old + x4) / 2.0 + 0.5 * h * x4d;

        if s < 0.09999 {
            s = 0.0;
            let xmy = disp * x3;
            let xic = x4 * disp / (vc * tgo);
            let tot = xmy + xic;
            let xmhe = -vm * he * x2;
            let xmhepz = xmhe / tp;

            array_tp.push(tp);
            array_tot.push(tot);
            array_xmhepz.push(xmhepz);
        }
    }

    Results {
        tp: array_tp,
        tot: array_tot,
        xmhepz: array_xmhepz,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data: Vec<Vec<f64>> = vec![
        results.tp.clone(),
        results.tot.clone(),
        results.xmhepz.clone(),
    ];
    let data_file = format!("{}/c20l7_datfil.txt", output_dir);
    save_data(&data_file, &data)?;

    println!("C20L7: Simulation complete");
    println!("  Output: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c20l7_runs() {
        let results = run();
        assert!(!results.tp.is_empty());
    }
}
