//! Chapter simulations from "Tactical and Strategic Missile Guidance, 7th ed"
//!
//! Each submodule contains lessons from a specific chapter.

pub mod c1;
pub mod c2;
pub mod c3;
pub mod c4;
pub mod c5;
pub mod c6;
pub mod c7;
pub mod c8;
pub mod c9;
pub mod c10;
pub mod c11;
pub mod c12;
pub mod c13;
pub mod c14;
pub mod c15;
pub mod c16;
pub mod c20;
pub mod c21;
pub mod c22;
pub mod c23;
pub mod c24;
pub mod c25;
pub mod c26;
pub mod c27;
pub mod c28;
pub mod c29;
pub mod c30;
pub mod c31;
pub mod c32;
pub mod c33;
pub mod c34;
pub mod c35;
pub mod c36;
pub mod c37;
pub mod c38;
pub mod c39;
pub mod c40;
pub mod c41;
pub mod c42;
pub mod c43;
pub mod c44;
pub mod c45;

// Placeholder modules (to be fully implemented)
pub mod c17;
pub mod c18;
pub mod c19;

/// List of all available simulations
pub fn list_simulations() -> Vec<(&'static str, &'static str)> {
    vec![
        ("c1l1", "Harmonic Oscillator - First Order Euler"),
        ("c1l2", "Harmonic Oscillator - Second Order RK2"),
        ("c1l3", "Digital Filter Step Response"),
        ("c2l1", "2D Tactical Missile-Target Engagement"),
        ("c2l2", "Linearized Engagement Model"),
        ("c3l1", "PN Miss Distance Analysis"),
        ("c4l1", "Gaussian Random Numbers"),
        ("c4l2", "Gaussian PDF Histogram"),
        ("c4l3", "Standard Deviation Calculation"),
        ("c4l4", "Low-Pass Filter White Noise"),
        ("c4l5", "Shaping Filter Monte Carlo"),
        ("c4l6", "Adjoint Shaping Filter"),
        ("c4l7", "Adjoint Noise Analysis"),
        ("c4l8", "Error Budget Analysis"),
        ("c5l1", "Fourth-Order Runge-Kutta"),
        ("c5l2", "Covariance Analysis"),
        ("c5l3", "Adjoint Method with Integration"),
        ("c6l1", "Higher-Order System Analysis"),
        ("c6l2", "Higher-Order with Integrators"),
        ("c6l3", "Fifth-Order Binomial Guidance"),
        ("c6l4", "Fifth-Order Adjoint Noise"),
        ("c7l1", "Fading Memory Filter"),
        ("c7l2", "Fading Memory Monte Carlo"),
        ("c7l3", "Fading Memory Adjoint"),
        ("c7l4", "Third-Order Digital Filter"),
        ("c8l1", "Augmented Proportional Navigation"),
        ("c8l2", "Optimal PN Adjoint"),
        ("c8l3", "OPN with TGO Error"),
        ("c9l1", "Kalman Filter Gains"),
        ("c9l2", "Kalman Filter Polynomial Model"),
        ("c9l3", "Kalman Filter Monte Carlo"),
        ("c9l4", "Variable Target Maneuver"),
        ("c9l5", "Optimal Guidance with Binomial Filter"),
        ("c10l1", "Ballistic Trajectory with Drag"),
        ("c11l1", "Reentry Target Trajectory"),
        ("c11l2", "2D Engagement with Reentry Target"),
        ("c12l1", "Extended Kalman Filter"),
        ("c12l2", "EKF Beta Estimation"),
        ("c12l3", "Polynomial Kalman Filter"),
        ("c13l1", "PN with Dynamics and Limits"),
        ("c14l1", "2D Engagement with Noise"),
        ("c14l2", "Command Guidance"),
        ("c15l1", "Flat Earth vs Spherical Earth"),
        ("c15l2", "Polar vs Cartesian Trajectory"),
        ("c15l3", "Circular Orbit"),
        ("c15l4", "Ballistic with Given Distance"),
        ("c15l5", "Minimum Energy Trajectory"),
        ("c15l6", "Flight Time vs Distance"),
        ("c15l7", "Velocity vs Flight Path Angle"),
        ("c15l8", "High Trajectory (Lob)"),
        ("c16l1", "Two-Stage Rocket Performance"),
        ("c16l2", "Two-Stage Rocket Trajectory"),
        ("c17l1", "Lambert Orbit Solver"),
        ("c17l2", "Lambert Direct Call"),
        ("c17l3", "Two-Stage with Lambert Guidance"),
        ("c17l4", "VTG-Based Guidance"),
        ("c17l5", "Ballistic Range Calculation"),
        ("c17l6", "Ballistic Propagation"),
        ("c18l1", "Strategic Intercept Ballistic"),
        ("c18l2", "Strategic Intercept Boosting"),
        ("c19l1", "Pulse Guidance"),
        ("c19l2", "Pulse Guidance Variant"),
        ("c20l1", "Target Jink/Displacement"),
        ("c20l2", "Target Displacement Adjoint"),
        ("c20l3", "Target Displacement Nonlinear"),
        ("c20l4", "Adjoint Third-Order Lag Filter"),
        ("c20l5", "Adjoint Second-Order Lag (Sweep TF)"),
        ("c20l6", "2D Engagement LOS Filter (Sweep THOM)"),
        ("c20l7", "Adjoint First-Order Lag (Heading Error)"),
        ("c20l8", "Adjoint Third-Order Lag (Total Miss)"),
        ("c20l9", "Adjoint Second-Order Lag (QSWITCH)"),
        ("c21l1", "Missile Aerodynamics"),
        ("c21l2", "Missile Aerodynamics - Transfer Function"),
        ("c22l1", "Missile Trim Analysis"),
        ("c22l2", "Rate Gyro Autopilot"),
        ("c22l3", "Frequency Response (Bode Plot)"),
        ("c22l4", "Rate Gyro with Actuator and Delay"),
        ("c22l5", "Numerical Frequency Response"),
        ("c23l1", "Three-Loop Autopilot"),
        ("c23l2", "Autopilot Frequency Response"),
        ("c23l3", "Radome Slope Analysis"),
        ("c23l4", "Nonlinear Autopilot Response"),
        ("c24l1", "Missile Aerodynamics (Flexible)"),
        ("c24l2", "Nonlinear Autopilot with Flexible Body"),
        ("c25l1", "Flexible Body Effects"),
        ("c25l2", "Flexible Body with Rate Feedback"),
        ("c25l3", "Flexible Body Two Modes with Notch"),
        ("c26l1", "Optimal Guidance - Riccati"),
        ("c26l2", "Optimal Control Simulation"),
        ("c26l3", "4th Order Riccati with Actuator"),
        ("c26l4", "Optimal Control with Gyro"),
        ("c26l5", "6th Order Riccati with Gyro"),
        ("c26l6", "Optimal Control Extended (6 gains)"),
        ("c26l7", "Frequency Response (Bode Plot)"),
        ("c26l8", "Frequency Response with Compensator"),
        ("c26l9", "Frequency Response with Gyro Comp"),
        ("c26l10", "Frequency Response Alternative"),
        ("c27l1", "Fading Memory Filter Miss"),
        ("c27l2", "Target Maneuver Miss with Fading Memory"),
        ("c27l3", "Miss for Various TGO Maneuver Start Times"),
        ("c27l4", "Target Maneuver Miss with Pre-computed Gains"),
        ("c27l5", "Miss for Various TGO with Blind Range"),
        ("c27l6", "Target Maneuver Miss with Gains and Blind Range"),
        ("c27l7", "Adjoint Model Target/Heading Error Miss"),
        ("c28l1", "Monte Carlo Miss Analysis"),
        ("c28l2", "Monte Carlo Miss with Shaping Filter"),
        ("c28l3", "Adjoint Model with Shaping Filter"),
        ("c29l1", "2D Weaving Target Engagement"),
        ("c29l2", "Target Weave Miss vs Flight Time"),
        ("c29l3", "Normalized Miss vs X Parameter"),
        ("c29l4", "Miss vs Flight Time (Guidance Laws)"),
        ("c29l5", "Adjoint Target Weave FFT Analysis"),
        ("c30l1", "Kalman Filter Acceleration Est"),
        ("c30l2", "Kalman Filter Singer Model"),
        ("c30l3", "Kalman Filter Frequency Est"),
        ("c31l1", "Multiple Model Adaptive Est"),
        ("c32l1", "Ballistic Target Intercept"),
        ("c32l2", "Ballistic Intercept Predictive"),
        ("c32l3", "Rolling Airframe (Open Loop)"),
        ("c32l4", "Rolling Airframe (Closed Loop)"),
        ("c33l1", "Predictive Guidance"),
        ("c34l1", "Shaping Filter Monte Carlo"),
        ("c34l2", "Adjoint Shaping Filter"),
        ("c34l3", "Kalman Filter Monte Carlo"),
        ("c35l1", "Optimal Guidance Gains"),
        ("c35l2", "Optimal Gains Oscillator"),
        ("c35l3", "Adjoint First-Order Lag"),
        ("c35l4", "Riccati with Autopilot"),
        ("c35l5", "Adjoint Autopilot Model"),
        ("c35l6", "Miss Distance Flight Times"),
        ("c36l1", "Impact Angle Control"),
        ("c36l2", "Impact Angle 2D Engagement"),
        ("c37l1", "Dive Guidance"),
        ("c37l2", "Polynomial Trajectory Matrix"),
        ("c37l3", "Dive Polynomial Guidance"),
        ("c38l1", "Trapezoidal Weave"),
        ("c38l2", "Monte Carlo Weave Target"),
        ("c38l3", "Adjoint Trapezoidal Weave"),
        ("c38l4", "Kalman Filter Weave Target"),
        ("c39l1", "Optimal Guidance with Autopilot"),
        ("c39l2", "Optimal Guidance Miss Analysis"),
        ("c39l3", "Optimal Guidance Full Autopilot"),
        ("c40l1", "3D Engagement"),
        ("c40l2", "3D Lambert Trajectory"),
        ("c40l3", "3D RK2 vs Kepler Comparison"),
        ("c40l4", "3D Strategic Intercept ZEM"),
        ("c40l5", "3D Reentry Guidance"),
        ("c41l1", "Kalman Filter Miss Est"),
        ("c41l2", "Kalman Filter Monte Carlo"),
        ("c42l1", "Optimal Guidance Monte Carlo"),
        ("c43l1", "Strategic Intercept"),
        ("c43l2", "Strategic Intercept Radar Track"),
        ("c44l1", "Ballistic Missile Defense"),
        ("c44l2", "BMD Flight Time vs Velocity"),
        ("c44l3", "BMD Interceptor Footprint"),
        ("c44l4", "BMD Target Location Coverage"),
        ("c44l5", "BMD Aim Point Coverage"),
        ("c45l1", "Theater Missile Defense"),
        ("c45l2", "TMD IRBM Trajectory"),
        ("c45l3", "TMD Intercept Simulation"),
        ("c45l4", "TMD Monte Carlo Simulation"),
    ]
}
