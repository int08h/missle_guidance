# Verification Summary

Last updated: January 30, 2026

## Overview

All 161 MATLAB chapter listings have been converted to Rust. Verification compares Rust output against MATLAB/Octave.

| Category           | Count   | Percentage   |
| ----------         | ------- | ------------ |
| PASS               |     102 |          63% |
| FAIL (Monte Carlo) |      22 |          14% |
| FAIL (Other)       |      21 |          13% |
| MATLAB No Output   |      16 |          10% |
| **Total**          | **161** | **100%**     |

**Deterministic Pass Rate: 70.3%** (102 of 145 comparable simulations)

## Passing Simulations (102)

All passing simulations have maximum relative error < 5% (most < 1e-6):

### Chapters 1-6
| Sim | Max Error | Description |
|-----|-----------|-------------|
| C1L1 | 4.45e-07 | Harmonic Oscillator - First Order Euler |
| C1L2 | 4.45e-07 | Harmonic Oscillator - Second Order RK2 |
| C1L3 | 5.02e-08 | Digital Filter Step Response |
| C2L1 | 4.43e-07 | 2D Tactical Missile-Target Engagement |
| C2L2 | 4.46e-07 | Linearized Engagement Model |
| C3L1 | 4.29e-07 | PN Miss Distance Analysis |
| C4L6 | 4.75e-07 | Adjoint Shaping Filter |
| C4L7 | 4.22e-07 | Adjoint Noise Analysis |
| C4L8 | 4.87e-07 | Error Budget Analysis |
| C5L1 | 3.63e-07 | Fourth-Order Runge-Kutta |
| C5L2 | 4.72e-07 | Covariance Analysis |
| C5L3 | 4.13e-07 | Adjoint Method with Integration |
| C6L1 | 4.21e-07 | Higher-Order System Analysis |
| C6L2 | 4.29e-07 | Higher-Order with Integrators |
| C6L3 | 3.30e-07 | Fifth-Order Binomial Guidance |
| C6L4 | 4.69e-07 | Fifth-Order Adjoint Noise |

### Chapters 8-16
| Sim | Max Error | Description |
|-----|-----------|-------------|
| C8L1 | 4.76e-07 | Augmented Proportional Navigation |
| C8L2 | 4.76e-07 | Optimal PN Adjoint |
| C8L3 | 4.27e-07 | PN with Noise |
| C9L1 | 3.85e-07 | Kalman Filter Gains |
| C9L5 | 4.13e-07 | Filter Analysis |
| C10L1 | 4.97e-07 | Ballistic Trajectory with Drag |
| C11L1 | 4.90e-07 | Reentry Target Trajectory |
| C13L1 | 4.31e-07 | PN with Dynamics and Limits |
| C14L2 | 4.47e-07 | Noise Analysis Adjoint |
| C15L1 | 4.45e-07 | Flat Earth vs Spherical Earth |
| C15L2 | 4.67e-07 | Weaving Target Analysis |
| C15L3 | 4.54e-07 | Multiple Weave Frequencies |
| C15L4 | 4.87e-07 | Extended Weave Analysis |
| C15L5 | 4.46e-03 | Nonlinear Weave |
| C15L6 | 1.21e-06 | Frequency Response |
| C15L7 | 4.72e-07 | Parameter Study |
| C15L8 | 1.88e-06 | Full Engagement |
| C16L1 | 4.70e-07 | Two-Stage Rocket Performance |
| C16L2 | 4.47e-07 | Radome Effects |

### Chapters 17-27
| Sim | Max Error | Description |
|-----|-----------|-------------|
| C17L5 | 1.54e-06 | 3D Geometry |
| C17L6 | 1.10e-05 | Full 3D Engagement |
| C20L1 | 4.84e-07 | Target Jink/Displacement |
| C20L2 | 4.96e-07 | Multiple Jink Patterns |
| C20L3 | 3.34e-03 | Parameter Sweep |
| C20L4 | 4.73e-07 | Noise Effects |
| C20L5 | 3.57e-07 | Filter Analysis |
| C20L6 | 2.22e-03 | Miss Statistics |
| C20L7 | 4.91e-07 | Extended Analysis |
| C20L8 | 4.70e-07 | Full Simulation |
| C20L9 | 4.83e-07 | Parameter Study |
| C21L1 | 4.67e-07 | Missile Aerodynamics |
| C21L2 | 4.89e-07 | Stability Analysis |
| C22L1 | 4.79e-07 | Missile Trim Analysis |
| C22L2 | 2.22e-07 | Lift Coefficient Study |
| C22L3 | 4.88e-07 | Drag Analysis |
| C22L4 | 4.70e-07 | Full Aero Model |
| C22L5 | 4.73e-07 | Parameter Sensitivity |
| C23L1 | 3.23e-07 | Three-Loop Autopilot |
| C23L2 | 4.01e-07 | Frequency Response |
| C23L3 | 4.58e-07 | Radome Slope |
| C23L4 | 4.82e-07 | Nonlinear Response |
| C24L1 | 4.67e-07 | Flexible Body Autopilot |
| C24L2 | 4.40e-07 | Flexible Body Nonlinear |
| C25L1 | 4.77e-07 | ZEM Guidance |
| C25L2 | 5.00e-07 | ZEM Extensions |
| C25L3 | 4.94e-07 | Full ZEM Analysis |
| C26L1 | 1.29e-03 | Optimal Guidance - Riccati |
| C26L2 | 4.93e-07 | Riccati Solutions |
| C26L4 | 4.90e-07 | Gain Schedules |
| C26L5 | 2.70e-03 | Time-Varying Gains |
| C26L6 | 4.93e-07 | Full Optimal |
| C26L7 | 4.74e-07 | Parameter Study |
| C26L8 | 4.06e-07 | Extended Analysis |
| C26L10 | 4.74e-07 | Comprehensive |
| C27L1 | 3.99e-07 | Fading Memory Filter |
| C27L4 | 3.80e-07 | Filter Tuning |
| C27L7 | 4.69e-07 | Full Filter Analysis |

### Chapters 28-45
| Sim | Max Error | Description |
|-----|-----------|-------------|
| C28L3 | 4.67e-07 | Orbit Mechanics |
| C29L1 | 3.72e-07 | 2D Weaving Target |
| C29L3 | 4.84e-07 | Weave Analysis |
| C29L4 | 4.04e-07 | Extended Weave |
| C32L1 | 4.53e-07 | Ballistic Intercept |
| C32L2 | 3.64e-07 | Intercept Analysis |
| C32L3 | 4.54e-07 | Full Trajectory |
| C32L4 | 4.89e-07 | Parameter Study |
| C33L1 | 4.53e-07 | Predictive Guidance |
| C34L2 | 3.90e-07 | Shaping Filter Adjoint |
| C35L1 | 4.90e-07 | Optimal Guidance Gains |
| C35L2 | 4.84e-07 | Gain Optimization |
| C35L3 | 4.76e-07 | Time-Optimal |
| C35L4 | 4.89e-07 | Energy-Optimal |
| C35L6 | 4.51e-07 | Comprehensive |
| C36L1 | 4.77e-07 | Impact Angle 2D |
| C36L2 | 4.63e-07 | Full Impact Angle |
| C37L1 | 4.74e-07 | Dive Guidance |
| C38L3 | 3.04e-03 | Trapezoidal Weave Variant |
| C39L1 | 4.92e-07 | Optimal with Autopilot |
| C39L2 | 4.22e-07 | Extended Optimal |
| C39L3 | 2.77e-04 | Full Analysis |
| C40L1 | 4.41e-07 | 3D Engagement |
| C41L1 | 2.03e-02 | Kalman Filter Estimation |
| C44L2 | 6.89e-08 | BMD Trajectory |
| C44L3 | 1.95e-07 | BMD Full |
| C44L4 | 1.20e-07 | BMD Extended |
| C44L5 | 1.94e-07 | BMD Comprehensive |
| C45L2 | 2.94e-03 | TMD Engagement |

## Monte Carlo Simulations (22 Expected Failures)

These simulations use `randn()` for random noise. Values differ each run:

| Sim | Notes |
|-----|-------|
| C4L1-L5 | Gaussian random number generation |
| C7L1, C7L3-L4 | Fading memory filter with noise |
| C12L1-L3 | Extended Kalman filter Monte Carlo |
| C14L1 | 2D engagement with noise |
| C28L1-L2 | Monte Carlo miss analysis |
| C30L1-L3 | Kalman acceleration estimation |
| C31L1 | Multiple model adaptive estimator |
| C34L1, C34L3 | Shaping filter Monte Carlo |
| C41L2 | Kalman filter Monte Carlo |
| C42L1 | Optimal guidance Monte Carlo |

## Other Failures (21)

These simulations have implementation differences or numerical issues:

| Sim | Max Error | Notes |
|-----|-----------|-------|
| C9L2-L4 | >1 | Filter implementation differences |
| C17L1 | 1.3 | 3D coordinate system |
| C27L2-L3, L5-L6 | >1 | Fading memory variants |
| C29L2, L5 | >1 | Weave target variants |
| C35L5 | 0.7 | Optimal gains variant |
| C38L1-L2, L4 | >0.2 | Trapezoidal weave |
| C40L4-L5 | 0.3 | 3D engagement variants |
| C43L1-L2 | >1 | Strategic intercept |
| C44L1 | 13 | BMD simplified output |
| C45L1, L3 | >1 | TMD variants |

## MATLAB Files Without Octave Output (16)

These files don't produce `datfil.txt` when run in Octave:

- C7L2 
- C11L2
- C17L2, C17L3, C17L4
- C18L1, C18L2
- C19L1, C19L2
- C26L3, C26L9
- C37L2, C37L3
- C40L2, C40L3
- C45L4

## Verification Method

1. Run Rust simulation: `./target/release/missile_guidance run <sim>`
2. Run MATLAB/Octave: `octave <SIM>.m`
3. Compare `datfil.txt` outputs using:
   - Relative error for values > 1e-6
   - Absolute error for near-zero values
   - Pass threshold: 5% maximum relative error

## Running Verification

```bash
./verify.sh
```

Results saved to `verification_results.txt`.

Note: verification of all files takes 15-20 minutes due to Monte Carlo simulations.
