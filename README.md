# Missile Guidance Simulations

Semi-complete Rust port of the MATLAB code from "Tactical and Strategic Missile Guidance, 7th ed" by Paul Zarchan.

## Available Simulations

161 simulations from 45 chapters:

- **C1-C3**: Numerical integration (Euler, RK2, RK4), basic engagement
- **C4-C6**: Monte Carlo analysis, adjoint methods, higher-order systems
- **C7-C9**: Fading memory filters, augmented PN, optimal guidance
- **C10-C14**: Lambert guidance, coasting flight, noise analysis
- **C15-C16**: Weaving targets, radome effects
- **C17-C19**: 3D engagements (not fully verified - Octave issues)
- **C20-C24**: Homing loops, aerodynamics, autopilot design
- **C25-C27**: ZEM guidance, optimal gains, filter design
- **C28-C31**: Monte Carlo miss analysis, ballistic prediction
- **C32-C34**: Extended Kalman filtering, shaping filters
- **C35-C39**: Optimal guidance gains, impact angle control
- **C40-C42**: 3D ballistic engagements, theater defense
- **C43-C45**: Strategic missile defense (BMD/TMD)

## Port/Verification Status

All 161 MATLAB simulations have been ported to Rust.

* 102 Rust ports pass verification and match the MATLAB/Octave output within a 5% relative error threshold.
* 22 Rust ports fail verification due to use of Monte Carlo and non-deterministic random numbers.
* 16 MATLAB originals produce no output. See `NO_OUTPUT.md` for details.
* 21 other verifications fail for various reasons. See `VERIFICATION_FAILS.md` for details.

| Status | Count | Description                        |
|--------|-------|------------------------------------|
| PASS             |     102 | Numerical results within tolerance |
| Monte Carlo FAIL |      22 | Expected - use random numbers      |
| MATLAB No Output |      16 | MATLAB issues, see `NO_OUTPUT.md`  |
| Other FAIL       |      21 | Various differences                |
| **Total**        | **161** |                                    |

### Buyer Beware

These are mechanical translations of MATLAB to Rust. There may be subtle differences in the numerical results. The
verification script compares the output of the Rust port to the original MATLAB output for assurance, but 100% 
accuracy is not guaranteed.

_Use at your own risk._

### Monte Carlo Simulations

These simulations use random numbers (`randn`) and produce different values each run:
- C4L1-L5, C7L1-L4, C12L1-L3, C14L1
- C28L1-L2, C30L1-L3, C31L1
- C34L1, C34L3, C41L2, C42L1

### MATLAB Files Without Octave Output

Some MATLAB files don't produce `datfil.txt` when run in Octave:
- C7L2, C11L2, C17L2-L4, C18L1-L2, C19L1-L2
- C26L3, C26L9, C37L2-L3, C40L2-L3, C45L4

See `NO_OUTPUT.md` for more details.

## Building

```bash
cargo build --release
```

## Usage

```bash
# List available simulations
./target/release/missile_guidance list

# Run a specific simulation
./target/release/missile_guidance run c2l1

# Run with custom output directory
./target/release/missile_guidance run c2l1 ./output

# Run all simulations
./target/release/missile_guidance run-all ./output
```

## Running Verification

```bash
# Full verification of all 161 simulations (requires Octave)
./verify.sh

# Verify a single simulation
./verify.sh c2l1     # Case insensitive
./verify.sh C10L1

# Results saved to verification_results.txt
```

## MATLAB Bug Fixes

The original MATLAB files had some bugs that were fixed:

- **C4L8.m**: Double plus `++` on lines 94-95. Fixed to use `+`.
- **C12L1.m**: Called `PROJECTC8L1` which doesn't exist. Fixed to call `PROJECTC12L1`.
- **C14L1.m**: Called `gaussc7` which doesn't exist. Fixed to use `SIGNOISE*randn`.
- **PROJECTC12L1.m**: Function name didn't match filename. Fixed function declaration.

## Output Format

Each simulation produces:
- `<sim>_datfil.txt` - ASCII data matching MATLAB format (space-separated scientific notation)
- `<sim>_*.png` - Plots (where applicable)

## Project Structure

```
src/
├── main.rs           # CLI entry point
├── lib.rs            # Library exports
├── plotting.rs       # PNG plot generation
├── chapters/         # Simulation implementations (c1/, c2/, ... c45/)
│   └── c{N}/l{M}.rs  # Chapter N, Listing M
└── utils/            # Shared utilities
    ├── kepler.rs     # Kepler orbit propagation
    ├── lambert3d.rs  # Lambert problem solver
    ├── predict.rs    # State prediction functions
    ├── project.rs    # Trajectory projection
    ├── gains.rs      # Guidance gain calculations
    └── rk2.rs        # Runge-Kutta integration
```

## License

Original MATLAB code from "Tactical and Strategic Missile Guidance, 7th ed" by Paul Zarchan. No licensing or copyright information is available.

Rust conversion copyright (c) 2026 Stuart Stock, all rights reserved. 

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
