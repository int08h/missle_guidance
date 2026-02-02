//! Chapter 20, Lesson 2: Target Displacement Adjoint
//!
//! Adjoint method for analyzing target displacement response.

use crate::save_data;

pub struct Results {
    pub tp: Vec<f64>,
    pub xmy: Vec<f64>,
    pub theory: Vec<f64>,
}

/// Run the C20L2 simulation
pub fn run() -> Results {
    let xnp: f64 = 3.0;
    let tau: f64 = 1.0;
    let tf: f64 = 5.0;
    let disp: f64 = 200.0;
    let vm: f64 = 3000.0;
    let _he = -disp / vm;
    let t: f64 = 0.0;
    let mut s: f64 = 0.0;
    let mut tp = t + 0.00001;
    let mut x1: f64 = 0.0;
    let mut x2: f64 = 0.0;
    let mut x3: f64 = 1.0;
    let mut x4: f64 = 0.0;
    let h: f64 = 0.01;

    let mut array_tp = Vec::new();
    let mut array_xmy = Vec::new();
    let mut array_theory = Vec::new();

    while tp <= tf - 0.00001 {
        s += h;
        let x1old = x1;
        let x2old = x2;
        let x3old = x3;
        let x4old = x4;

        // First Euler step
        let mut x1d = x2;
        let mut x2d = x3;
        let y1 = (x4 - x2) / tau;
        let tgo = tp + 0.00001;
        let mut x3d = xnp * y1 / tgo;
        let mut x4d = -y1;

        x1 += h * x1d;
        x2 += h * x2d;
        x3 += h * x3d;
        x4 += h * x4d;
        tp += h;

        // Second evaluation for RK2
        x1d = x2;
        x2d = x3;
        let y1 = (x4 - x2) / tau;
        let tgo = tp + 0.00001;
        x3d = xnp * y1 / tgo;
        x4d = -y1;

        // RK2 averaging
        x1 = (x1old + x1) / 2.0 + 0.5 * h * x1d;
        x2 = (x2old + x2) / 2.0 + 0.5 * h * x2d;
        x3 = (x3old + x3) / 2.0 + 0.5 * h * x3d;
        x4 = (x4old + x4) / 2.0 + 0.5 * h * x4d;

        if s < 0.09999 {
            s = 0.0;
            let xmy = disp * x3;
            let x = tgo / tau;
            let theory = disp * (-x).exp() * (1.0 - 2.0 * x + 0.5 * x * x);

            array_tp.push(tp);
            array_xmy.push(xmy);
            array_theory.push(theory);
        }
    }

    Results {
        tp: array_tp,
        xmy: array_xmy,
        theory: array_theory,
    }
}

pub fn run_and_save(output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    let data: Vec<Vec<f64>> = vec![
        results.tp.clone(),
        results.xmy.clone(),
        results.theory.clone(),
    ];
    let data_file = format!("{}/c20l2_datfil.txt", output_dir);
    save_data(&data_file, &data)?;

    println!("C20L2: Simulation complete");
    println!("  Output: {}", data_file);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c20l2_runs() {
        let results = run();
        assert!(!results.tp.is_empty());
    }
}
