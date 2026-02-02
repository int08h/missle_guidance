# Why MATLAB Files Don't Produce Output in Octave

This document summarizes the investigation into why certain MATLAB files don't produce `datfil.txt` output when run in Octave.

## Categories of Issues Found

| Category | Count | Files |
|----------|-------|-------|
| **Graphics blocking** | 12 | C11L2, C17L3, C24L1, C26L1, C26L4, C26L6, C27L6, C35L3, C35L4, C36L1, C45L4 |
| **Missing save command** | 3 | C17L2, C37L2, C40L3 |
| **Undefined function/variable** | 6 | C7L2, C18L1, C18L2, C26L3, C26L9, C40L2 |
| **Syntax errors/bugs** | 3 | C26L2, C26L5, C37L3 |
| **Octave syntax incompatibility** | 1 | C29L2 |
| **Loop/initialization issues** | 6 | C17L4, C19L1, C19L2, C35L1, C35L5, C39L1 |
| **Performance/timeout** | 1 | C35L2 |

## Detailed Breakdown

### 1. Graphics Commands Blocking (Most Common)

Scripts have `figure`/`plot` commands before the `save` statement. In headless Octave, these fail or hang, preventing the save from executing.

| File | Issue |
|------|-------|
| C11L2 | `figure`/`plot` commands at lines 133-140 before save at line 143 |
| C17L3 | Graphics operations at lines 136-137 before save at line 142 |
| C24L1 | Graphics commands at lines 82-89 before save at line 92 |
| C26L1 | Multiple `figure`/`plot` commands at lines 58-65 before save at line 55 |
| C26L4 | Graphics commands at lines 73-86 before save at line 89 |
| C26L6 | Three `figure`/`plot` blocks at lines 67-79 before save at line 83 |
| C27L6 | `figure`/`plot` at lines 207-208 before save at line 214 |
| C35L3 | Graphics commands at lines 62-71 after save but `clc` may interfere |
| C35L4 | Graphics commands at lines 82-86 after save at line 79 |
| C36L1 | Multiple graphics commands at lines 62-78 before save at line 81 |
| C45L4 | `figure`/`plot` at lines 519-524 before save at line 530 |

### 2. Missing Save Command

Scripts only display results to console, never write to files.

| File | Purpose |
|------|---------|
| C17L2 | Lambert solver demonstration - displays velocities to console only (20 lines) |
| C37L2 | Polynomial coefficient calculator - console output only, no file I/O |
| C40L3 | Kepler orbit verification script - displays comparison results to console |

### 3. Undefined Functions/Variables (MATLAB Bugs)

| File | Bug | Fix |
|------|-----|-----|
| C7L2 | Line 58: Calls `gaussc7(SIGNOISE)` which doesn't exist | Replace with `SIGNOISE*randn` |
| C18L1 | Line 46: Calls `lambertpz()` but function file defines `lambert()` | Fix function name |
| C18L2 | Lines 62, 65: Calls `predictb()` and `lambertpz()` which may not exist | Verify function availability |
| C26L3 | Line 41: Uses `GT` but only `G` is defined | Change to `G'` (transpose) |
| C26L9 | Line 28: References `ArrayPHASE` which is never computed | Add phase calculation |
| C40L2 | Line 78: Calls `distance3d()` but function is `distance3dkm()` | Fix function name |

### 4. Syntax Errors/Bugs

| File | Line | Bug | Fix |
|------|------|-----|-----|
| C26L2 | 47 | `if S>=.0-9999` evaluates to `S>=-9999` (always true) | Should be `if S>=0.01` |
| C26L5 | 8 | `KDEL=9000.;;` has double semicolon | Remove extra semicolon |
| C37L3 | 44 | `while RM@ >= 0` has invalid `@` symbol | Should be `while RM1 >= 0` |

### 5. Octave Syntax Incompatibility

| File | Line | Issue | Fix |
|------|------|-------|-----|
| C29L2 | 66 | Uses `save datfil.txt output /ascii` (MATLAB syntax) | Change `/ascii` to `-ascii` |

### 6. Loop/Initialization Issues

| File | Issue |
|------|-------|
| C17L4 | Depends on `lambertpz` function from C17L2 which itself has issues |
| C19L1 | Loop condition `while VC>=0` likely false initially, loop never executes |
| C19L2 | Arrays only initialized inside loop; if conditions not met, arrays undefined |
| C35L1 | Dynamic array allocation without pre-allocation may fail in strict Octave |
| C35L5 | Numerical instability: division by near-zero (`D=X+GAM` where GAM=0.00001) |
| C39L1 | Floating-point threshold `S1>=.009999` may have precision issues |

### 7. Performance/Timeout

| File | Issue |
|------|-------|
| C35L2 | ~100,000 iterations of 4x4 matrix operations; likely times out in interpreted Octave |

## Recommendations

### Quick Fixes (Could Enable Verification)

1. **C29L2**: Change `/ascii` to `-ascii` on line 66
2. **C26L2**: Fix line 47 from `if S>=.0-9999` to `if S>=0.01`
3. **C26L5**: Remove extra semicolon on line 8
4. **C37L3**: Fix line 44 from `RM@` to `RM1`
5. **C26L3**: Change line 41 from `GT` to `G'`
6. **C7L2**: Replace `gaussc7(SIGNOISE)` with `SIGNOISE*randn` on line 58

### Requires More Work

- Graphics-blocking files need graphics commands moved after save or wrapped in conditionals
- Missing save commands need explicit file output added
- Function name mismatches need coordination with utility files

### Cannot Be Fixed (By Design)

- Console-only demonstration scripts (C17L2, C40L3) - intentionally interactive
- Performance-limited scripts (C35L2) - need MATLAB or compiled code
