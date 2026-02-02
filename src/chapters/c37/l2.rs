//! Chapter 37, Lesson 2: Polynomial Trajectory Matrix Calculation
//!
//! Computes X2 coefficient for polynomial guidance using matrix inversion.
//! This is a simple matrix calculation, not a time-domain simulation.
//!
//! NOTE: The MATLAB code does NOT produce datfil output - it only prints values
//! to the console. This is a computation-only listing.

pub struct Results {
    pub x2sim3: f64,
    pub x2theory3: f64,
    pub x2sim4: f64,
    pub x2theory4: f64,
    pub x2sim5: f64,
    pub x2theory5: f64,
}

/// Simple matrix inversion for small matrices using Gaussian elimination
#[allow(clippy::needless_range_loop)]
fn invert_matrix(a: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = a.len();
    let mut aug = vec![vec![0.0; 2 * n]; n];

    // Create augmented matrix [A | I]
    for i in 0..n {
        for j in 0..n {
            aug[i][j] = a[i][j];
        }
        aug[i][n + i] = 1.0;
    }

    // Gaussian elimination with partial pivoting
    for i in 0..n {
        // Find pivot
        let mut max_row = i;
        for k in i + 1..n {
            if aug[k][i].abs() > aug[max_row][i].abs() {
                max_row = k;
            }
        }
        aug.swap(i, max_row);

        // Scale row i
        let pivot = aug[i][i];
        if pivot.abs() < 1e-15 {
            continue;
        }
        for j in 0..2 * n {
            aug[i][j] /= pivot;
        }

        // Eliminate column i
        for k in 0..n {
            if k != i {
                let factor = aug[k][i];
                for j in 0..2 * n {
                    aug[k][j] -= factor * aug[i][j];
                }
            }
        }
    }

    // Extract inverse matrix
    let mut inv = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            inv[i][j] = aug[i][n + j];
        }
    }
    inv
}

/// Matrix-vector multiplication
fn mat_vec_mul(a: &[Vec<f64>], x: &[f64]) -> Vec<f64> {
    let n = a.len();
    let mut result = vec![0.0; n];
    for i in 0..n {
        for j in 0..a[i].len() {
            result[i] += a[i][j] * x[j];
        }
    }
    result
}

/// Run the C37L2 simulation
pub fn run() -> Results {
    let r0: f64 = 10000.0;
    let vm: f64 = 250.0;
    let tfdes: f64 = 50.0;
    let e0: f64 = 0.0;

    // N=3 case (4x4 matrix)
    let a3 = vec![
        vec![1.0, 0.0, 0.0, 0.0],
        vec![0.0, 1.0, 0.0, 0.0],
        vec![1.0, tfdes, tfdes.powi(2), tfdes.powi(3)],
        vec![0.0, 1.0, 2.0 * tfdes, 3.0 * tfdes.powi(2)],
    ];
    let x0_3 = vec![r0, -vm * e0.cos(), 0.0, -vm];
    let a3_inv = invert_matrix(&a3);
    let xpz1 = mat_vec_mul(&a3_inv, &x0_3);
    let x2sim3 = xpz1[2];
    let n3: f64 = 3.0;
    let x2theory3 = (n3 - 1.0) * ((n3 - 2.0 + 2.0 * e0.cos()) * vm * tfdes - n3 * r0)
        / (2.0 * tfdes * tfdes);

    // N=4 case (5x5 matrix)
    let a4 = vec![
        vec![1.0, 0.0, 0.0, 0.0, 0.0],
        vec![0.0, 1.0, 0.0, 0.0, 0.0],
        vec![1.0, tfdes, tfdes.powi(2), tfdes.powi(3), tfdes.powi(4)],
        vec![0.0, 1.0, 2.0 * tfdes, 3.0 * tfdes.powi(2), 4.0 * tfdes.powi(3)],
        vec![0.0, 0.0, 2.0, 6.0 * tfdes, 12.0 * tfdes.powi(2)],
    ];
    let x0_4 = vec![r0, -vm * e0.cos(), 0.0, -vm, 0.0];
    let a4_inv = invert_matrix(&a4);
    let xpz2 = mat_vec_mul(&a4_inv, &x0_4);
    let x2sim4 = xpz2[2];
    let n4: f64 = 4.0;
    let x2theory4 = (n4 - 1.0) * ((n4 - 2.0 + 2.0 * e0.cos()) * vm * tfdes - n4 * r0)
        / (2.0 * tfdes * tfdes);

    // N=5 case (6x6 matrix)
    let a5 = vec![
        vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        vec![0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
        vec![1.0, tfdes, tfdes.powi(2), tfdes.powi(3), tfdes.powi(4), tfdes.powi(5)],
        vec![0.0, 1.0, 2.0 * tfdes, 3.0 * tfdes.powi(2), 4.0 * tfdes.powi(3), 5.0 * tfdes.powi(4)],
        vec![0.0, 0.0, 2.0, 6.0 * tfdes, 12.0 * tfdes.powi(2), 20.0 * tfdes.powi(3)],
        vec![0.0, 0.0, 0.0, 6.0, 24.0 * tfdes, 60.0 * tfdes.powi(2)],
    ];
    let x0_5 = vec![r0, -vm * e0.cos(), 0.0, -vm, 0.0, 0.0];
    let a5_inv = invert_matrix(&a5);
    let xpz3 = mat_vec_mul(&a5_inv, &x0_5);
    let x2sim5 = xpz3[2];
    let n5: f64 = 5.0;
    let x2theory5 = (n5 - 1.0) * ((n5 - 2.0 + 2.0 * e0.cos()) * vm * tfdes - n5 * r0)
        / (2.0 * tfdes * tfdes);

    Results {
        x2sim3,
        x2theory3,
        x2sim4,
        x2theory4,
        x2sim5,
        x2theory5,
    }
}

pub fn run_and_save(_output_dir: &str) -> std::io::Result<Results> {
    let results = run();

    // NOTE: MATLAB's C37L2 does NOT produce datfil.txt output.
    // It only prints values to the console. We match this behavior.
    println!("C37L2: Simulation finished");
    println!("  X2SIM3 = {:.6e}, X2THEORY3 = {:.6e}", results.x2sim3, results.x2theory3);
    println!("  X2SIM4 = {:.6e}, X2THEORY4 = {:.6e}", results.x2sim4, results.x2theory4);
    println!("  X2SIM5 = {:.6e}, X2THEORY5 = {:.6e}", results.x2sim5, results.x2theory5);

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c37l2_runs() {
        let results = run();
        // Check that simulation matches theory
        assert!((results.x2sim3 - results.x2theory3).abs() < 1e-6);
        assert!((results.x2sim4 - results.x2theory4).abs() < 1e-6);
        assert!((results.x2sim5 - results.x2theory5).abs() < 1e-6);
    }
}
