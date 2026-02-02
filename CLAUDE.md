# CLAUDE.md

## Project Overview

Rust implementation of missile guidance simulations from "Tactical and Strategic Missile Guidance, 7th ed" by Paul Zarchan. The goal is to produce numerically identical output to the original MATLAB code.

**Status: All 161 MATLAB chapter listings have been converted to Rust.**

## Build and Run

```bash
cargo build --release
./target/release/missile_guidance run c2l1        # Run single simulation
./target/release/missile_guidance run-all ./output # Run all simulations
./target/release/missile_guidance list            # List available
```

## Verification

```bash
./verify.sh          # Compare all Rust outputs against MATLAB/Octave (requires Octave)
./verify.sh c2l1     # Verify a single simulation (case insensitive)
./verify.sh C10L1    # Verify a single simulation
```

## Project Structure

```
src/
├── main.rs           # CLI entry point
├── lib.rs            # Library exports
├── plotting.rs       # PNG plot generation
├── chapters/         # Simulation implementations (c1/, c2/, ... c45/)
│   └── c{N}/l{M}.rs  # Chapter N, Listing M (161 files total)
└── utils/            # Shared utilities
    ├── kepler.rs     # Kepler orbit propagation
    ├── lambert3d.rs  # Lambert problem solver
    ├── predict.rs    # State prediction functions
    ├── project.rs    # Trajectory projection
    ├── gains.rs      # Guidance gain calculations
    └── rk2.rs        # Runge-Kutta integration
```

## Conversion Status

| Metric | Count |
|--------|-------|
| Total MATLAB listings | 161 |
| Rust implementations | 161 |
| **Conversion complete** | **100%** |

## Verification Results

| Category | Count | Notes |
|----------|-------|-------|
| PASS | 102 | Numerical match within 5% tolerance |
| Monte Carlo FAIL | 22 | Expected - use random numbers |
| Other FAIL | 21 | Implementation differences |
| MATLAB No Output | 16 | Octave compatibility issues |

**Deterministic Pass Rate: 70.3%** (102 of 145 comparable simulations)

## Output Matching Requirements

**Rust output should match MATLAB output exactly.**

### Acceptable Differences

- Relative error < 5% for deterministic simulations
- Monte Carlo simulations: same structure (rows/columns), different values due to `randn`
- Near-zero values: absolute error < 1e-10

### Monte Carlo Simulations

These use random numbers and will differ each run:
- C4L1-L5, C7L1-L4, C9L2-L4, C12L1-L3, C14L1
- C28L1-L2, C30L1-L3, C31L1
- C34L1, C34L3, C38L2, C38L4, C41L2, C42L1, C43L2

### Common Mismatch Causes

- **Integration step timing**: MATLAB uses specific step sizes (H) and sampling periods (TS)
- **TGO updates**: Time-to-go must update at each RK substep
- **Output columns**: Must output same columns in same order as MATLAB
- **Floating-point order**: Operations must follow MATLAB's evaluation order

## MATLAB Bugs Found and Fixed

These bugs were in the original MATLAB files:

- **C4L8.m**: Double plus `++` on lines 94-95. Fixed to use `+`.
- **C12L1.m**: Called `PROJECTC8L1` (doesn't exist). Fixed to `PROJECTC12L1`.
- **C14L1.m**: Called `gaussc7` (doesn't exist). Fixed to `SIGNOISE*randn`.
- **PROJECTC12L1.m**: Function name didn't match filename. Fixed declaration.

## MATLAB Files Without Octave Output

Some MATLAB files don't produce `datfil.txt` when run in Octave:
- C7L2, C11L2, C17L2-L4, C18L1-L2, C19L1-L2
- C26L3, C26L9, C37L2-L3, C40L2-L3, C45L4

These likely require MATLAB-specific features or have graphics-only output.

## Utility Functions

MATLAB utility functions have been implemented in Rust:

| Rust File | MATLAB Functions |
|-----------|------------------|
| `kepler.rs` | kepler.m, KEPLER1.m |
| `lambert3d.rs` | LAMBERT3D.m |
| `lambert2d.rs` | lamberpz.m |
| `predict.rs` | PREDICT1.m, predict44.m, predict45.m, predictb.m, predictg.m, predictpz.m, distance3dkm.m, LAUNCHLOGIC.m |
| `project.rs` | PROJECT.m, PROJECT34.m, PROJECT6S.m, PROJECTC12L1.m, PROJECTC12L2.m |
| `gains.rs` | GENERATEGAINS.m |
| `initial.rs` | INITIALPZ.m |

## Code Quality

- Clippy clean (no warnings)
- 197 unit tests passing
- All simulations produce output

See `VERIFICATION_SUMMARY.md` for verification results and `NO_OUTPUT.md` for MATLAB files without output.