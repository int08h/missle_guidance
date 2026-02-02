//! State projection functions for dynamics simulation
//!
//! Implements PROJECT, PROJECT34, PROJECT6S and chapter-specific variants

/// Result from PROJECT function
#[derive(Debug, Clone, Copy)]
pub struct ProjectResult {
    pub yb: f64,
    pub ydb: f64,
    pub ytddb: f64,
    pub ytdddb: f64,
}

/// Project state forward with 4th-order dynamics (PROJECT.m)
///
/// Projects position, velocity, acceleration, and jerk forward in time
/// with a simple oscillator model.
///
/// # Arguments
/// * `tp` - Current time (unused in calculation)
/// * `ts` - Time step to project
/// * `yph` - Initial position
/// * `ydph` - Initial velocity
/// * `ytddph` - Initial acceleration
/// * `ytdddph` - Initial jerk
/// * `hp` - Integration step size
/// * `xnlp` - Nonlinear acceleration offset
/// * `wph` - Natural frequency
#[allow(clippy::too_many_arguments)]
pub fn project(
    _tp: f64,
    ts: f64,
    yph: f64,
    ydph: f64,
    ytddph: f64,
    ytdddph: f64,
    hp: f64,
    xnlp: f64,
    wph: f64,
) -> ProjectResult {
    let mut t = 0.0;
    let mut y = yph;
    let mut yd = ydph;
    let mut ytdd = ytddph;
    let mut ytddd = ytdddph;
    let w = wph;
    let xnl = xnlp;
    let h = hp;

    while t <= ts - 0.0001 {
        let ytdddd = -w * w * ytdd;
        ytddd += h * ytdddd;
        ytdd += h * ytddd;
        let ydd = ytdd - xnl;
        yd += h * ydd;
        y += h * yd;
        t += h;
    }

    ProjectResult {
        yb: y,
        ydb: yd,
        ytddb: ytdd,
        ytdddb: ytddd,
    }
}

/// Result from PROJECT34 function
#[derive(Debug, Clone, Copy)]
pub struct Project34Result {
    pub y1: f64,
    pub y2: f64,
    pub y3: f64,
}

/// Project 3-state system with optional Poisson dynamics (PROJECT34.m)
///
/// # Arguments
/// * `ts` - Time step to project
/// * `y1_init` - Initial state 1
/// * `y2_init` - Initial state 2
/// * `y3_init` - Initial state 3
/// * `hp` - Integration step size
/// * `xnl` - Nonlinear offset
/// * `poession` - Poisson parameter (for exponential decay)
pub fn project34(
    ts: f64,
    y1_init: f64,
    y2_init: f64,
    y3_init: f64,
    hp: f64,
    xnl: f64,
    poession: f64,
) -> Project34Result {
    let mut t = 0.0;
    let mut y1 = y1_init;
    let mut y2 = y2_init;
    let mut y3 = y3_init;
    let h = hp;

    while t <= ts - 0.0001 {
        let y3d = -poession * y3;
        y3 += h * y3d;
        let y2d = y3 - xnl;
        y2 += h * y2d;
        let y1d = y2;
        y1 += h * y1d;
        t += h;
    }

    Project34Result { y1, y2, y3 }
}

/// Result from PROJECT6S function
#[derive(Debug, Clone, Copy)]
pub struct Project6SResult {
    pub y: [f64; 6],
}

/// Project 6-state system with coupled oscillatory dynamics (PROJECT6S.m)
///
/// # Arguments
/// * `ts` - Time step to project
/// * `y_init` - Initial 6-state vector
/// * `hp` - Integration step size
/// * `w` - Natural frequency
/// * `zeta` - Damping ratio
pub fn project6s(
    ts: f64,
    y_init: &[f64; 6],
    hp: f64,
    w: f64,
    zeta: f64,
) -> Project6SResult {
    let mut t = 0.0;
    let mut y = *y_init;
    let h = hp;

    while t <= ts - 0.0001 {
        // 6th order dynamics with oscillatory coupling
        let y6d = -2.0 * zeta * w * y[5] - w * w * y[4];
        y[5] += h * y6d;
        y[4] += h * y[5];
        let y4d = y[4];
        y[3] += h * y4d;
        let y3d = y[3];
        y[2] += h * y3d;
        let y2d = y[2];
        y[1] += h * y2d;
        let y1d = y[1];
        y[0] += h * y1d;
        t += h;
    }

    Project6SResult { y }
}

/// Chapter 12 Lesson 1 specific projection function
pub fn project_c12l1(
    ts: f64,
    y_init: f64,
    yd_init: f64,
    ydd_init: f64,
    hp: f64,
    xnc: f64,
    tau: f64,
) -> (f64, f64, f64) {
    let mut t = 0.0;
    let mut y = y_init;
    let mut yd = yd_init;
    let mut ydd = ydd_init;
    let h = hp;

    while t <= ts - 0.0001 {
        let yddd = (xnc - ydd) / tau;
        ydd += h * yddd;
        yd += h * ydd;
        y += h * yd;
        t += h;
    }

    (y, yd, ydd)
}

/// Chapter 12 Lesson 2 specific projection function with damping
#[allow(clippy::too_many_arguments)]
pub fn project_c12l2(
    ts: f64,
    y_init: f64,
    yd_init: f64,
    ydd_init: f64,
    yddd_init: f64,
    hp: f64,
    xnc: f64,
    w: f64,
    zeta: f64,
) -> (f64, f64, f64, f64) {
    let mut t = 0.0;
    let mut y = y_init;
    let mut yd = yd_init;
    let mut ydd = ydd_init;
    let mut yddd = yddd_init;
    let h = hp;

    while t <= ts - 0.0001 {
        let ydddd = w * w * (xnc - ydd) - 2.0 * zeta * w * yddd;
        yddd += h * ydddd;
        ydd += h * yddd;
        yd += h * ydd;
        y += h * yd;
        t += h;
    }

    (y, yd, ydd, yddd)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_project_basic() {
        let result = project(0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.01, 0.0, 1.0);
        // With zero initial acceleration and unity velocity, position should increase
        assert!(result.yb > 0.0);
    }

    #[test]
    fn test_project34_basic() {
        let result = project34(1.0, 0.0, 1.0, 0.0, 0.01, 0.0, 0.0);
        assert!(result.y1 > 0.0);
    }
}
