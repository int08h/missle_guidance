//! Kalman filter gain generation
//!
//! Implements GENERATEGAINS for computing time-varying Kalman filter gains

/// Container for generated gains over time
#[derive(Debug, Clone)]
pub struct GeneratedGains {
    pub c1: Vec<f64>,
    pub c2: Vec<f64>,
    pub c3: Vec<f64>,
    pub c4: Vec<f64>,
    pub c5: Vec<f64>,
    pub c6: Vec<f64>,
}

/// Generate time-varying Kalman filter gains (GENERATEGAINS.m)
///
/// Solves the Riccati equation backwards in time to generate optimal
/// filter gains for state estimation.
///
/// # Arguments
/// * `tau` - Time constant
/// * `w` - Natural frequency
/// * `z` - Damping ratio
/// * `wz` - Secondary frequency
/// * `gam` - Measurement noise factor
/// * `tf` - Final time
/// * `ts` - Sampling interval for output
///
/// # Returns
/// Generated gains at each sampling instant
pub fn generate_gains(
    tau: f64,
    w: f64,
    z: f64,
    wz: f64,
    gam: f64,
    tf: f64,
    ts: f64,
) -> GeneratedGains {
    // Initialize matrices
    // F matrix (6x6) - system dynamics
    let mut f = [[0.0f64; 6]; 6];
    f[0][1] = 1.0;
    f[1][2] = 1.0;
    f[1][3] = -1.0;
    f[1][5] = 1.0 / (wz * wz);
    f[3][4] = 1.0;
    f[4][5] = 1.0;
    f[5][3] = -w * w / tau;
    f[5][4] = -w * w * (2.0 * z / w + tau) / tau;
    f[5][5] = -w * w * (1.0 / (w * w) + 2.0 * z * tau / w) / tau;

    // G vector (6x1) - input
    let mut g = [0.0f64; 6];
    g[5] = w * w / tau;

    // S matrix (6x6) - Riccati solution, initialized with measurement
    let mut s = [[0.0f64; 6]; 6];
    s[0][0] = 1.0;

    let h = 0.0001;
    let mut t = 0.0;
    let mut s1 = 0.0;

    let mut c1 = Vec::new();
    let mut c2 = Vec::new();
    let mut c3 = Vec::new();
    let mut c4 = Vec::new();
    let mut c5 = Vec::new();
    let mut c6 = Vec::new();

    while t < tf - 0.0001 {
        s1 += h;
        let sold = s;

        // Compute S*F
        let mut sf = [[0.0f64; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                for k in 0..6 {
                    sf[i][j] += s[i][k] * f[k][j];
                }
            }
        }

        // Compute G' * S
        let mut gts = [0.0f64; 6];
        for j in 0..6 {
            for i in 0..6 {
                gts[j] += g[i] * s[i][j];
            }
        }

        // Compute C = (1/gamma) * G' * S
        let gam_inv = 1.0 / gam;
        let c: Vec<f64> = gts.iter().map(|&x| gam_inv * x).collect();

        // Compute C' * C
        let mut ctc = [[0.0f64; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                ctc[i][j] = c[i] * c[j];
            }
        }

        // Compute gamma * C' * C
        let mut ctbc = [[0.0f64; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                ctbc[i][j] = gam * ctc[i][j];
            }
        }

        // Compute SD = S*F + F'*S - gamma*C'*C
        let mut sd = [[0.0f64; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                sd[i][j] = sf[i][j] + sf[j][i] - ctbc[i][j];
            }
        }

        // Euler step
        for i in 0..6 {
            for j in 0..6 {
                s[i][j] += h * sd[i][j];
            }
        }
        t += h;

        // RK2 correction
        // Recompute sf, gts, c, ctc, ctbc, sd with new s
        let mut sf2 = [[0.0f64; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                for k in 0..6 {
                    sf2[i][j] += s[i][k] * f[k][j];
                }
            }
        }

        let mut gts2 = [0.0f64; 6];
        for j in 0..6 {
            for i in 0..6 {
                gts2[j] += g[i] * s[i][j];
            }
        }

        let c2_temp: Vec<f64> = gts2.iter().map(|&x| gam_inv * x).collect();

        let mut ctc2 = [[0.0f64; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                ctc2[i][j] = c2_temp[i] * c2_temp[j];
            }
        }

        let mut ctbc2 = [[0.0f64; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                ctbc2[i][j] = gam * ctc2[i][j];
            }
        }

        let mut sd2 = [[0.0f64; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                sd2[i][j] = sf2[i][j] + sf2[j][i] - ctbc2[i][j];
            }
        }

        // RK2 averaging
        for i in 0..6 {
            for j in 0..6 {
                s[i][j] = 0.5 * (sold[i][j] + s[i][j]) + 0.5 * h * sd2[i][j];
            }
        }

        // Store gains at sampling intervals
        if s1 >= ts - 0.0001 {
            s1 = 0.0;

            // Recompute c for storage
            let mut gts_final = [0.0f64; 6];
            for j in 0..6 {
                for i in 0..6 {
                    gts_final[j] += g[i] * s[i][j];
                }
            }
            let c_final: Vec<f64> = gts_final.iter().map(|&x| -gam_inv * x).collect();

            c1.push(c_final[0]);
            c2.push(c_final[1]);
            c3.push(c_final[2]);
            c4.push(c_final[3]);
            c5.push(c_final[4]);
            c6.push(c_final[5]);
        }
    }

    GeneratedGains { c1, c2, c3, c4, c5, c6 }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generate_gains() {
        let gains = generate_gains(1.0, 10.0, 0.7, 5.0, 1.0, 10.0, 0.1);
        assert!(!gains.c1.is_empty());
        assert_eq!(gains.c1.len(), gains.c2.len());
    }
}
