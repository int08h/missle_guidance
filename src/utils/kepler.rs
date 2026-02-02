//! Kepler orbit propagation solvers
//!
//! Contains both the simple kepler solver (kepler.m) and the advanced
//! universal variable formulation solver (KEPLER1.m)

use super::constants::{GM_FT, GM_KM, EARTH_RADIUS_KM};

/// Result from simple Kepler solver
#[derive(Debug, Clone, Copy)]
pub struct KeplerResult {
    pub xf: f64,
    pub yf: f64,
    pub rf: f64,
}

/// Kepler solver using binary search (kepler.m)
///
/// Solves for the final position given initial conditions and desired flight time.
///
/// # Arguments
/// * `gam_deg` - Flight path angle in degrees
/// * `v` - Initial velocity
/// * `tf_des` - Desired time of flight
/// * `xlong_m_deg` - Initial longitude in degrees
/// * `alt_m` - Initial altitude
///
/// # Returns
/// `KeplerResult` containing final position and radius
pub fn kepler_solve(
    gam_deg: f64,
    v: f64,
    tf_des: f64,
    xlong_m_deg: f64,
    alt_m: f64,
) -> KeplerResult {
    let pi = std::f64::consts::PI;
    let deg_rad = 360.0 / (2.0 * pi);
    let a = 2.0926e7; // Earth radius in ft
    let gm = GM_FT;

    let gam = gam_deg / deg_rad;
    let xlong_m = xlong_m_deg / deg_rad;
    let xm = (a + alt_m) * xlong_m.cos();
    let ym = (a + alt_m) * xlong_m.sin();
    let r0 = (xm * xm + ym * ym).sqrt();

    let mut icount = 0;
    let mut phi_max = pi / 2.0;
    let mut phi_min = 0.0;
    let mut phi = 45.0 / deg_rad;
    let mut tf = 100000.0;
    let mut phi_old = phi;
    let mut told = tf;
    let mut rf = 0.0;

    while (tf_des - tf).abs() > 0.00000001 * tf_des {
        let _phi_deg = phi * deg_rad;
        let cphi = phi.cos();
        let sphi = phi.sin();

        let top = v * v * r0 * r0 * gam.cos() * gam.cos();
        let bot = gm * (1.0 - cphi) + r0 * v * v * gam.cos() * (phi + gam).cos();
        rf = top / bot;

        let xlam = r0 * v * v / gm;
        let top1 = gam.tan() * (1.0 - cphi) + (1.0 - xlam) * sphi;
        let bot1p = (1.0 - cphi) / (xlam * gam.cos() * gam.cos());
        let bot1 = (2.0 - xlam) * (bot1p + (gam + phi).cos() / gam.cos());
        let top2 = 2.0 * gam.cos();
        let inner = 2.0 / xlam - 1.0;

        if inner <= 0.0 {
            // Adjust bounds if hyperbolic
            phi_max = phi;
            phi = (phi_max + phi_min) / 2.0;
            continue;
        }

        let bot2 = xlam * inner.powf(1.5);
        let top3 = inner.sqrt();
        let bot3 = gam.cos() / (phi / 2.0).tan() - gam.sin();
        let temp = (top2 / bot2) * top3.atan2(bot3);
        tf = r0 * (top1 / bot1 + temp) / (v * gam.cos());

        icount += 1;

        if tf > tf_des {
            phi_max = phi;
        } else {
            phi_min = phi;
        }

        let xnext = if icount == 1 {
            (phi_max + phi_min) / 2.0
        } else {
            let mut xnext = phi + (phi - phi_old) * (tf_des - tf) / (tf - told);
            if xnext > phi_max || xnext < phi_min {
                xnext = (phi_max + phi_min) / 2.0;
            }
            xnext
        };

        phi_old = phi;
        told = tf;
        phi = xnext;

        if icount > 100 {
            break;
        }
    }

    let xf = rf * phi.cos();
    let yf = rf * phi.sin();

    KeplerResult { xf, yf, rf }
}

/// State vector for Kepler1 solver
#[derive(Debug, Clone, Copy)]
pub struct StateVector {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub vx: f64,
    pub vy: f64,
    pub vz: f64,
}

impl StateVector {
    pub fn new(x: f64, y: f64, z: f64, vx: f64, vy: f64, vz: f64) -> Self {
        Self { x, y, z, vx, vy, vz }
    }

    pub fn from_array(arr: &[f64; 6]) -> Self {
        Self {
            x: arr[0], y: arr[1], z: arr[2],
            vx: arr[3], vy: arr[4], vz: arr[5],
        }
    }

    pub fn to_array(self) -> [f64; 6] {
        [self.x, self.y, self.z, self.vx, self.vy, self.vz]
    }
}

/// Advanced Kepler solver using universal variable formulation (KEPLER1.m)
///
/// Propagates a state vector forward in time using the universal Kepler equation.
/// Works for all orbit types (elliptic, parabolic, hyperbolic).
///
/// # Arguments
/// * `x0` - Initial state vector [x, y, z, vx, vy, vz] in km and km/s
/// * `t0` - Initial time in seconds
/// * `t1` - Final time in seconds
///
/// # Returns
/// Final state vector at time t1
pub fn kepler1(x0: &StateVector, t0: f64, t1: f64) -> StateVector {
    let gmx = GM_KM;
    let rex = EARTH_RADIUS_KM;
    let tlimit = 1.0e-10;
    let kn = 10;
    let muqr = 1.0;

    let dt_orig = t1 - t0;

    if dt_orig.abs() < tlimit {
        return *x0;
    }

    let time_factor = (rex.powi(3) / gmx).sqrt();
    let vex = rex / time_factor;
    let dt = dt_orig / time_factor;

    // Normalize state
    let r0 = [x0.x / rex, x0.y / rex, x0.z / rex];
    let v0 = [x0.vx / vex, x0.vy / vex, x0.vz / vex];

    let r0_mag = (r0[0].powi(2) + r0[1].powi(2) + r0[2].powi(2)).sqrt();
    let v0_mag = (v0[0].powi(2) + v0[1].powi(2) + v0[2].powi(2)).sqrt();
    let d0 = r0[0] * v0[0] + r0[1] * v0[1] + r0[2] * v0[2];
    let sigma0 = d0 / muqr;
    let alp0 = 2.0 / r0_mag - v0_mag * v0_mag;

    let _a0 = if alp0 == 0.0 { 1.0e30 } else { 1.0 / alp0 };

    let mut x = if alp0 <= 0.0 {
        0.1 * dt / r0_mag
    } else {
        alp0 * dt
    };

    let mut cy;
    let mut sy;
    let mut u1 = 0.0;
    let mut u2 = 0.0;
    let mut u3 = 0.0;
    let mut dfx = r0_mag;

    for _k in 0..kn {
        let y = alp0 * x * x;

        if alp0 < 0.0 {
            let yqr = (-y).sqrt();
            cy = (1.0 - yqr.cosh()) / y;
            sy = (yqr.sinh() - yqr) / yqr.powi(3);
        } else if alp0 == 0.0 {
            cy = 0.5;
            sy = 1.0 / 6.0;
        } else {
            let yqr = y.sqrt();
            cy = (1.0 - yqr.cos()) / y;
            sy = (yqr - yqr.sin()) / yqr.powi(3);
        }

        u1 = x * (1.0 - y * sy);
        u2 = x * x * cy;
        u3 = x * x * x * sy;

        let fx = r0_mag * u1 + sigma0 * u2 + u3 - dt * muqr;
        dfx = sigma0 * u1 + (1.0 - alp0 * r0_mag) * u2 + r0_mag;
        let dfx2 = sigma0 * (1.0 - y * cy) + (1.0 - alp0 * r0_mag) * u1;

        let sdfx = dfx.signum();
        let dx2 = 16.0 * dfx * dfx - 20.0 * fx * dfx2;

        let dx = if dx2 > 0.0 {
            5.0 * fx / (dfx + sdfx * dx2.sqrt())
        } else {
            0.5 * x
        };

        x -= dx;
    }

    let rmag = dfx;
    let f = 1.0 - u2 / r0_mag;
    let g = dt - u3 / muqr;
    let df = -muqr * u1 / (rmag * r0_mag);
    let dg = 1.0 - u2 / rmag;

    StateVector {
        x: (f * r0[0] + g * v0[0]) * rex,
        y: (f * r0[1] + g * v0[1]) * rex,
        z: (f * r0[2] + g * v0[2]) * rex,
        vx: (df * r0[0] + dg * v0[0]) * vex,
        vy: (df * r0[1] + dg * v0[1]) * vex,
        vz: (df * r0[2] + dg * v0[2]) * vex,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kepler_solve() {
        let result = kepler_solve(30.0, 24000.0, 2200.0, 0.0, 0.0);
        assert!(result.rf > 0.0);
        assert!(result.xf != 0.0 || result.yf != 0.0);
    }

    #[test]
    fn test_kepler1_circular_orbit() {
        // Test with a circular orbit
        let rex = EARTH_RADIUS_KM;
        let v_circ = (GM_KM / rex).sqrt();

        let x0 = StateVector::new(rex, 0.0, 0.0, 0.0, v_circ, 0.0);
        let x1 = kepler1(&x0, 0.0, 1000.0);

        // Check conservation of radius (approximately)
        let r1 = (x1.x.powi(2) + x1.y.powi(2) + x1.z.powi(2)).sqrt();
        assert!((r1 - rex).abs() / rex < 0.01);
    }
}
