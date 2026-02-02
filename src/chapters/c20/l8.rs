//! Chapter 20, Listing 8: Adjoint with Third-Order Lag (Total Miss)
//!
//! Simulates adjoint response with third-order lag filter (8 states).
//! Computes total miss as sum of displacement and ZEM terms.

use crate::save_data;

pub struct Results {
    pub tp: Vec<f64>,
    pub tot: Vec<f64>,
}

/// Run the C20L8 simulation
pub fn run() -> Results {
    let xnp: f64 = 3.0;
    let tau: f64 = 1.0;
    let tf: f64 = 10.0;
    let disp: f64 = 1.0;
    let vc: f64 = 4000.0;
    let vm: f64 = 3000.0;
    let _he: f64 = -disp / vm;

    let mut _t: f64 = 0.0;
    let mut s: f64 = 0.0;
    let mut tp: f64 = _t + 0.00001;
    let mut x2: f64 = 0.0;
    let mut x3: f64 = 1.0;
    let mut x4: f64 = 0.0;
    let mut x5: f64 = 0.0;
    let mut x6: f64 = 0.0;
    let mut x7: f64 = 0.0;
    let mut x8: f64 = 0.0;
    let h: f64 = 0.01;

    let mut array_tp = Vec::new();
    let mut array_tot = Vec::new();

    while tp <= tf - 0.00001 {
        s += h;
        let x2old = x2;
        let x3old = x3;
        let x4old = x4;
        let x5old = x5;
        let x6old = x6;
        let x7old = x7;
        let x8old = x8;

        // First derivative evaluation
        let mut x2d = x3;
        let mut y1 = 5.0 * (5.0 * x5 / tau + x4) / tau;
        let mut tgo = tp + 0.00001;
        let mut x3d = y1 / (vc * tgo);
        let mut x4d = -y1;
        let mut x5d = -5.0 * x5 / tau + 5.0 * x6 * xnp * vc / tau;
        let mut x6d = -5.0 * x6 / tau + 5.0 * x7 / tau;
        let mut x7d = -5.0 * x7 / tau + 5.0 * x8 / tau;
        let mut x8d = -5.0 * x8 / tau - x2;

        // Euler step
        x2 += h * x2d;
        x3 += h * x3d;
        x4 += h * x4d;
        x5 += h * x5d;
        x6 += h * x6d;
        x7 += h * x7d;
        x8 += h * x8d;
        tp += h;

        // Second derivative evaluation
        x2d = x3;
        y1 = 5.0 * (5.0 * x5 / tau + x4) / tau;
        tgo = tp + 0.00001;
        x3d = y1 / (vc * tgo);
        x4d = -y1;
        x5d = -5.0 * x5 / tau + 5.0 * x6 * xnp * vc / tau;
        x6d = -5.0 * x6 / tau + 5.0 * x7 / tau;
        x7d = -5.0 * x7 / tau + 5.0 * x8 / tau;
        x8d = -5.0 * x8 / tau - x2;

        // RK2 averaging
        x2 = (x2old + x2) / 2.0 + 0.5 * h * x2d;
        x3 = (x3old + x3) / 2.0 + 0.5 * h * x3d;
        x4 = (x4old + x4) / 2.0 + 0.5 * h * x4d;
        x5 = (x5old + x5) / 2.0 + 0.5 * h * x5d;
        x6 = (x6old + x6) / 2.0 + 0.5 * h * x6d;
        x7 = (x7old + x7) / 2.0 + 0.5 * h * x7d;
        x8 = (x8old + x8) / 2.0 + 0.5 * h * x8d;

        if s < 0.09999 {
            s = 0.0;
            let xmy = disp * x3;
            let xic = x4 * disp / (vc * tgo);
            let tot = xmy + xic;

            array_tp.push(tp);
            array_tot.push(tot);
        }
    }

    Results {
        tp: array_tp,
        tot: array_tot,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data: Vec<Vec<f64>> = vec![
        results.tp.clone(),
        results.tot.clone(),
    ];
    let data_file = format!("{}/c20l8_datfil.txt", output_dir);
    save_data(&data_file, &data)?;

    println!("C20L8: Simulation complete");
    println!("  Output: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c20l8_runs() {
        let results = run();
        assert!(!results.tp.is_empty());
    }
}
