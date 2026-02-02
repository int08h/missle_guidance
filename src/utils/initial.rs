//! Initial condition calculators
//!
//! Implements INITIALPZ for computing ballistic trajectory initial conditions

use super::constants::G_ACCEL;

/// Result from initial trajectory calculation
#[derive(Debug, Clone, Copy)]
pub struct InitialResult {
    pub rt1f: f64,
    pub rt2f: f64,
    pub tfdes: f64,
}

/// Calculate initial trajectory conditions with drag (INITIALPZ.m)
///
/// Propagates trajectory forward until a specified altitude is reached,
/// accounting for atmospheric drag.
///
/// # Arguments
/// * `rt2des` - Desired final altitude
/// * `rt1ic` - Initial downrange position
/// * `rt2ic` - Initial altitude
/// * `vt1ic` - Initial downrange velocity
/// * `vt2ic` - Initial vertical velocity
/// * `beta` - Ballistic coefficient
///
/// # Returns
/// Final position and time to reach desired altitude
pub fn initial_pz(
    rt2des: f64,
    rt1ic: f64,
    rt2ic: f64,
    vt1ic: f64,
    vt2ic: f64,
    beta: f64,
) -> InitialResult {
    let mut rt1 = rt1ic;
    let mut rt2 = rt2ic;
    let mut vt1 = vt1ic;
    let mut vt2 = vt2ic;
    let mut t = 0.0;
    let h = 0.01;

    while rt2 > rt2des {
        let rt1old = rt1;
        let rt2old = rt2;
        let vt1old = vt1;
        let vt2old = vt2;

        // Compute atmospheric density
        let rho = if rt2 <= 30000.0 {
            0.002378 * (-rt2 / 30000.0).exp()
        } else {
            0.0034 * (-rt2 / 22000.0).exp()
        };

        // Velocity magnitude and dynamic pressure
        let vt = (vt1 * vt1 + vt2 * vt2).sqrt();
        let q = 0.5 * rho * vt * vt;

        // Flight path angle
        let gamt = (-vt2).atan2(vt1);

        // Acceleration with drag
        let at1 = -G_ACCEL * q * gamt.cos() / beta;
        let at2 = -G_ACCEL + G_ACCEL * q * gamt.sin() / beta;

        // Euler step
        rt1 += h * vt1;
        rt2 += h * vt2;
        vt1 += h * at1;
        vt2 += h * at2;
        t += h;

        // Second evaluation for RK2
        let rho = if rt2 <= 30000.0 {
            0.002378 * (-rt2 / 30000.0).exp()
        } else {
            0.0034 * (-rt2 / 22000.0).exp()
        };
        let vt = (vt1 * vt1 + vt2 * vt2).sqrt();
        let q = 0.5 * rho * vt * vt;
        let gamt = (-vt2).atan2(vt1);
        let at1 = -G_ACCEL * q * gamt.cos() / beta;
        let at2 = -G_ACCEL + G_ACCEL * q * gamt.sin() / beta;

        // RK2 averaging
        rt1 = 0.5 * (rt1old + rt1 + h * vt1);
        rt2 = 0.5 * (rt2old + rt2 + h * vt2);
        vt1 = 0.5 * (vt1old + vt1 + h * at1);
        vt2 = 0.5 * (vt2old + vt2 + h * at2);
    }

    InitialResult {
        rt1f: rt1,
        rt2f: rt2,
        tfdes: t,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_initial_pz() {
        // Simple falling object
        let result = initial_pz(0.0, 0.0, 10000.0, 1000.0, 0.0, 1000.0);
        assert!(result.tfdes > 0.0);
        assert!(result.rt1f > 0.0); // Should have moved downrange
    }
}
