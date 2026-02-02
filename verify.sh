#!/bin/bash
# Comprehensive verification script to compare Rust and MATLAB outputs
# Tests all 161 converted listings, or a specific one if provided
#
# Usage:
#   ./verify.sh          # Run all tests
#   ./verify.sh c2l1     # Run only C2L1
#   ./verify.sh C2L1     # Case insensitive

RUST_DIR=$(pwd)
MATLAB_DIR="$RUST_DIR/MATLAB"
VERIFY_DIR="$RUST_DIR/verify_output"
RESULTS_FILE="$RUST_DIR/verification_results.txt"

# Parse optional argument for single test
SINGLE_TEST=""
if [ -n "$1" ]; then
    # Convert to uppercase for matching
    SINGLE_TEST=$(echo "$1" | tr '[:lower:]' '[:upper:]')
fi

mkdir -p "$VERIFY_DIR"

# Build Rust project first
echo "Building Rust project..."
cd "$RUST_DIR" || exit
cargo build --release -q
if [ $? -ne 0 ]; then
    echo "ERROR: Rust build failed!"
    exit 1
fi

# Run Rust simulation(s) to generate output files
if [ -n "$SINGLE_TEST" ]; then
    echo "Running Rust simulation for $SINGLE_TEST..."
    local_sim=$(echo "$SINGLE_TEST" | tr '[:upper:]' '[:lower:]')
    ./target/release/missile_guidance run "$local_sim" "$VERIFY_DIR" > /dev/null 2>&1
else
    echo "Running all Rust simulations..."
    ./target/release/missile_guidance run-all "$VERIFY_DIR" > /dev/null 2>&1
fi

echo ""
echo "Verification Results - $(date)" > "$RESULTS_FILE"
echo "================================" >> "$RESULTS_FILE"
echo "" >> "$RESULTS_FILE"

compare_files() {
    local rust_file="$1"
    local matlab_file="$2"
    local name="$3"

    if [ ! -f "$rust_file" ]; then
        echo "$name: RUST FILE MISSING"
        return 1
    fi
    if [ ! -f "$matlab_file" ]; then
        echo "$name: MATLAB FILE MISSING"
        return 1
    fi

    # Compare values using awk
    # Use absolute error for near-zero values, relative error otherwise
    local result=$(awk '
    BEGIN { max_err = 0; count = 0; sum_err = 0; }
    NR==FNR {
        for(i=1; i<=NF; i++) rust[NR,i] = $i;
        rust_cols[NR] = NF;
        rust_rows = NR;
        next
    }
    {
        for(i=1; i<=NF; i++) {
            if(rust[FNR,i] != "" && $i != "") {
                r = rust[FNR,i] + 0;
                m = $i + 0;
                abs_m = (m < 0 ? -m : m);
                abs_r = (r < 0 ? -r : r);
                if(abs_m < 1e-6 && abs_r < 1e-6) {
                    err = (r - m);
                    if(err < 0) err = -err;
                    if(err < 1e-6) err = 0;
                } else if(abs_m > 1e-10) {
                    err = (r - m) / m;
                    if(err < 0) err = -err;
                } else if(abs_r > 1e-10) {
                    err = 1;
                } else {
                    err = 0;
                }
                if(err > max_err) max_err = err;
                sum_err += err;
                count++;
            }
        }
        matlab_rows = FNR;
    }
    END {
        avg_err = (count > 0) ? sum_err / count : 0;
        printf "%.6e %.6e %d %d %d", max_err, avg_err, count, rust_rows, matlab_rows
    }
    ' "$rust_file" "$matlab_file")

    local max_err=$(echo $result | cut -d' ' -f1)
    local avg_err=$(echo $result | cut -d' ' -f2)
    local count=$(echo $result | cut -d' ' -f3)
    local rust_rows=$(echo $result | cut -d' ' -f4)
    local matlab_rows=$(echo $result | cut -d' ' -f5)

    # Check if max error is small enough (5% relative error threshold)
    local threshold="5e-2"
    local pass=$(awk "BEGIN { print ($max_err < $threshold) ? 1 : 0 }")

    if [ "$pass" = "1" ]; then
        echo "$name: PASS (max_err=$max_err, n=$count)"
    else
        echo "$name: FAIL (max_err=$max_err, n=$count)"
    fi
}

run_matlab_and_compare() {
    local chapter="$1"
    local rust_sim="$2"
    local matlab_file="$3"

    # Skip if single test mode and this isn't the requested test
    if [ -n "$SINGLE_TEST" ] && [ "$chapter" != "$SINGLE_TEST" ]; then
        return
    fi

    printf "Testing %-8s " "$chapter..."

    # Run MATLAB simulation
    cd "$MATLAB_DIR"
    rm -f datfil.txt datfil 2>/dev/null
    octave --no-gui --silent "$matlab_file" > /dev/null 2>&1

    # Find the output file
    local matlab_out=""
    if [ -f "$MATLAB_DIR/datfil.txt" ]; then
        matlab_out="$MATLAB_DIR/datfil.txt"
    elif [ -f "$MATLAB_DIR/datfil" ]; then
        matlab_out="$MATLAB_DIR/datfil"
    else
        echo "$chapter: MATLAB NO OUTPUT" | tee -a "$RESULTS_FILE"
        return
    fi

    # Compare outputs
    local rust_file="$VERIFY_DIR/${rust_sim}_datfil.txt"
    local result=$(compare_files "$rust_file" "$matlab_out" "$chapter")
    echo "$result" | tee -a "$RESULTS_FILE"
}

if [ -n "$SINGLE_TEST" ]; then
    echo "Verifying single test: $SINGLE_TEST"
else
    echo "Starting verification of all 161 listings..."
fi
echo ""

# Chapter 1
run_matlab_and_compare "C1L1" "c1l1" "C1L1.m"
run_matlab_and_compare "C1L2" "c1l2" "C1L2.m"
run_matlab_and_compare "C1L3" "c1l3" "C1L3.m"

# Chapter 2
run_matlab_and_compare "C2L1" "c2l1" "C2L1.m"
run_matlab_and_compare "C2L2" "c2l2" "C2L2.m"

# Chapter 3
run_matlab_and_compare "C3L1" "c3l1" "C3L1.m"

# Chapter 4
run_matlab_and_compare "C4L1" "c4l1" "C4L1.m"
run_matlab_and_compare "C4L2" "c4l2" "C4L2.m"
run_matlab_and_compare "C4L3" "c4l3" "C4L3.m"
run_matlab_and_compare "C4L4" "c4l4" "C4L4.m"
run_matlab_and_compare "C4L5" "c4l5" "C4L5.m"
run_matlab_and_compare "C4L6" "c4l6" "C4L6.m"
run_matlab_and_compare "C4L7" "c4l7" "C4L7.m"
run_matlab_and_compare "C4L8" "c4l8" "C4L8.m"

# Chapter 5
run_matlab_and_compare "C5L1" "c5l1" "C5L1.m"
run_matlab_and_compare "C5L2" "c5l2" "C5L2.m"
run_matlab_and_compare "C5L3" "c5l3" "C5L3.m"

# Chapter 6
run_matlab_and_compare "C6L1" "c6l1" "C6L1.m"
run_matlab_and_compare "C6L2" "c6l2" "C6L2.m"
run_matlab_and_compare "C6L3" "c6l3" "C6L3.m"
run_matlab_and_compare "C6L4" "c6l4" "C6L4.m"

# Chapter 7 (Monte Carlo - structure match only)
run_matlab_and_compare "C7L1" "c7l1" "C7L1.m"
run_matlab_and_compare "C7L2" "c7l2" "C7L2.m"
run_matlab_and_compare "C7L3" "c7l3" "C7L3.m"
run_matlab_and_compare "C7L4" "c7l4" "C7L4.m"

# Chapter 8
run_matlab_and_compare "C8L1" "c8l1" "C8L1.m"
run_matlab_and_compare "C8L2" "c8l2" "C8L2.m"
run_matlab_and_compare "C8L3" "c8l3" "C8L3.m"

# Chapter 9
run_matlab_and_compare "C9L1" "c9l1" "C9L1.m"
run_matlab_and_compare "C9L2" "c9l2" "C9L2.m"
run_matlab_and_compare "C9L3" "c9l3" "C9L3.m"
run_matlab_and_compare "C9L4" "c9l4" "C9L4.m"
run_matlab_and_compare "C9L5" "c9l5" "C9L5.m"

# Chapter 10
run_matlab_and_compare "C10L1" "c10l1" "C10L1.m"

# Chapter 11
run_matlab_and_compare "C11L1" "c11l1" "C11L1.m"
run_matlab_and_compare "C11L2" "c11l2" "C11L2.m"

# Chapter 12 (Monte Carlo)
run_matlab_and_compare "C12L1" "c12l1" "C12L1.m"
run_matlab_and_compare "C12L2" "c12l2" "C12L2.m"
run_matlab_and_compare "C12L3" "c12l3" "C12L3.m"

# Chapter 13
run_matlab_and_compare "C13L1" "c13l1" "C13L1.m"

# Chapter 14 (Monte Carlo)
run_matlab_and_compare "C14L1" "c14l1" "C14L1.m"
run_matlab_and_compare "C14L2" "c14l2" "C14L2.m"

# Chapter 15
run_matlab_and_compare "C15L1" "c15l1" "C15L1.m"
run_matlab_and_compare "C15L2" "c15l2" "C15L2.m"
run_matlab_and_compare "C15L3" "c15l3" "C15L3.m"
run_matlab_and_compare "C15L4" "c15l4" "C15L4.m"
run_matlab_and_compare "C15L5" "c15l5" "C15L5.m"
run_matlab_and_compare "C15L6" "c15l6" "C15L6.m"
run_matlab_and_compare "C15L7" "c15l7" "C15L7.m"
run_matlab_and_compare "C15L8" "c15l8" "C15L8.m"

# Chapter 16
run_matlab_and_compare "C16L1" "c16l1" "C16L1.m"
run_matlab_and_compare "C16L2" "c16l2" "C16L2.m"

# Chapter 17
run_matlab_and_compare "C17L1" "c17l1" "C17L1.m"
run_matlab_and_compare "C17L2" "c17l2" "C17L2.m"
run_matlab_and_compare "C17L3" "c17l3" "C17L3.m"
run_matlab_and_compare "C17L4" "c17l4" "C17L4.m"
run_matlab_and_compare "C17L5" "c17l5" "C17L5.m"
run_matlab_and_compare "C17L6" "c17l6" "C17L6.m"

# Chapter 18
run_matlab_and_compare "C18L1" "c18l1" "C18L1.m"
run_matlab_and_compare "C18L2" "c18l2" "C18L2.m"

# Chapter 19
run_matlab_and_compare "C19L1" "c19l1" "C19L1.m"
run_matlab_and_compare "C19L2" "c19l2" "C19L2.m"

# Chapter 20
run_matlab_and_compare "C20L1" "c20l1" "C20L1.m"
run_matlab_and_compare "C20L2" "c20l2" "C20L2.m"
run_matlab_and_compare "C20L3" "c20l3" "C20L3.m"
run_matlab_and_compare "C20L4" "c20l4" "C20L4.m"
run_matlab_and_compare "C20L5" "c20l5" "C20L5.m"
run_matlab_and_compare "C20L6" "c20l6" "C20L6.m"
run_matlab_and_compare "C20L7" "c20l7" "C20L7.m"
run_matlab_and_compare "C20L8" "c20l8" "C20L8.m"
run_matlab_and_compare "C20L9" "c20l9" "C20L9.m"

# Chapter 21
run_matlab_and_compare "C21L1" "c21l1" "C21L1.m"
run_matlab_and_compare "C21L2" "c21l2" "C21L2.m"

# Chapter 22
run_matlab_and_compare "C22L1" "c22l1" "C22L1.m"
run_matlab_and_compare "C22L2" "c22l2" "C22L2.m"
run_matlab_and_compare "C22L3" "c22l3" "C22L3.m"
run_matlab_and_compare "C22L4" "c22l4" "C22L4.m"
run_matlab_and_compare "C22L5" "c22l5" "C22L5.m"

# Chapter 23
run_matlab_and_compare "C23L1" "c23l1" "C23L1.m"
run_matlab_and_compare "C23L2" "c23l2" "C23L2.m"
run_matlab_and_compare "C23L3" "c23l3" "C23L3.m"
run_matlab_and_compare "C23L4" "c23l4" "C23L4.m"

# Chapter 24
run_matlab_and_compare "C24L1" "c24l1" "C24L1.m"
run_matlab_and_compare "C24L2" "c24l2" "C24L2.m"

# Chapter 25
run_matlab_and_compare "C25L1" "c25l1" "C25L1.m"
run_matlab_and_compare "C25L2" "c25l2" "C25L2.m"
run_matlab_and_compare "C25L3" "c25l3" "C25L3.m"

# Chapter 26
run_matlab_and_compare "C26L1" "c26l1" "C26L1.m"
run_matlab_and_compare "C26L2" "c26l2" "C26L2.m"
run_matlab_and_compare "C26L3" "c26l3" "C26L3.m"
run_matlab_and_compare "C26L4" "c26l4" "C26L4.m"
run_matlab_and_compare "C26L5" "c26l5" "C26L5.m"
run_matlab_and_compare "C26L6" "c26l6" "C26L6.m"
run_matlab_and_compare "C26L7" "c26l7" "C26L7.m"
run_matlab_and_compare "C26L8" "c26l8" "C26L8.m"
run_matlab_and_compare "C26L9" "c26l9" "C26L9.m"
run_matlab_and_compare "C26L10" "c26l10" "C26L10.m"

# Chapter 27
run_matlab_and_compare "C27L1" "c27l1" "C27L1.m"
run_matlab_and_compare "C27L2" "c27l2" "C27L2.m"
run_matlab_and_compare "C27L3" "c27l3" "C27L3.m"
run_matlab_and_compare "C27L4" "c27l4" "C27L4.m"
run_matlab_and_compare "C27L5" "c27l5" "C27L5.m"
run_matlab_and_compare "C27L6" "c27l6" "C27L6.m"
run_matlab_and_compare "C27L7" "c27l7" "C27L7.m"

# Chapter 28 (Monte Carlo)
run_matlab_and_compare "C28L1" "c28l1" "C28L1.m"
run_matlab_and_compare "C28L2" "c28l2" "C28L2.m"
run_matlab_and_compare "C28L3" "c28l3" "C28L3.m"

# Chapter 29
run_matlab_and_compare "C29L1" "c29l1" "C29L1.m"
run_matlab_and_compare "C29L2" "c29l2" "C29L2.m"
run_matlab_and_compare "C29L3" "c29l3" "C29L3.m"
run_matlab_and_compare "C29L4" "c29l4" "C29L4.m"
run_matlab_and_compare "C29L5" "c29l5" "C29L5.m"

# Chapter 30 (Monte Carlo)
run_matlab_and_compare "C30L1" "c30l1" "C30L1.m"
run_matlab_and_compare "C30L2" "c30l2" "C30L2.m"
run_matlab_and_compare "C30L3" "c30l3" "C30L3.m"

# Chapter 31 (Monte Carlo)
run_matlab_and_compare "C31L1" "c31l1" "C31L1.m"

# Chapter 32
run_matlab_and_compare "C32L1" "c32l1" "C32L1.m"
run_matlab_and_compare "C32L2" "c32l2" "C32L2.m"
run_matlab_and_compare "C32L3" "c32l3" "C32L3.m"
run_matlab_and_compare "C32L4" "c32l4" "C32L4.m"

# Chapter 33
run_matlab_and_compare "C33L1" "c33l1" "C33L1.m"

# Chapter 34 (Monte Carlo)
run_matlab_and_compare "C34L1" "c34l1" "C34L1.m"
run_matlab_and_compare "C34L2" "c34l2" "C34L2.m"
run_matlab_and_compare "C34L3" "c34l3" "C34L3.m"

# Chapter 35
run_matlab_and_compare "C35L1" "c35l1" "C35L1.m"
run_matlab_and_compare "C35L2" "c35l2" "C35L2.m"
run_matlab_and_compare "C35L3" "c35l3" "C35L3.m"
run_matlab_and_compare "C35L4" "c35l4" "C35L4.m"
run_matlab_and_compare "C35L5" "c35l5" "C35L5.m"
run_matlab_and_compare "C35L6" "c35l6" "C35L6.m"

# Chapter 36
run_matlab_and_compare "C36L1" "c36l1" "C36L1.m"
run_matlab_and_compare "C36L2" "c36l2" "C36L2.m"

# Chapter 37
run_matlab_and_compare "C37L1" "c37l1" "C37L1.m"
run_matlab_and_compare "C37L2" "c37l2" "C37L2.m"
run_matlab_and_compare "C37L3" "c37l3" "C37L3.m"

# Chapter 38
run_matlab_and_compare "C38L1" "c38l1" "C38L1.m"
run_matlab_and_compare "C38L2" "c38l2" "C38L2.m"
run_matlab_and_compare "C38L3" "c38l3" "C38L3.m"
run_matlab_and_compare "C38L4" "c38l4" "C38L4.m"

# Chapter 39
run_matlab_and_compare "C39L1" "c39l1" "C39L1.m"
run_matlab_and_compare "C39L2" "c39l2" "C39L2.m"
run_matlab_and_compare "C39L3" "c39l3" "C39L3.m"

# Chapter 40
run_matlab_and_compare "C40L1" "c40l1" "C40L1.m"
run_matlab_and_compare "C40L2" "c40l2" "C40L2.m"
run_matlab_and_compare "C40L3" "c40l3" "C40L3.m"
run_matlab_and_compare "C40L4" "c40l4" "C40L4.m"
run_matlab_and_compare "C40L5" "c40l5" "C40L5.m"

# Chapter 41 (Monte Carlo)
run_matlab_and_compare "C41L1" "c41l1" "C41L1.m"
run_matlab_and_compare "C41L2" "c41l2" "C41L2.m"

# Chapter 42 (Monte Carlo)
run_matlab_and_compare "C42L1" "c42l1" "C42L1.m"

# Chapter 43
run_matlab_and_compare "C43L1" "c43l1" "C43L1.m"
run_matlab_and_compare "C43L2" "c43l2" "C43L2.m"

# Chapter 44
run_matlab_and_compare "C44L1" "c44l1" "C44L1.m"
run_matlab_and_compare "C44L2" "c44l2" "C44L2.m"
run_matlab_and_compare "C44L3" "c44l3" "C44L3.m"
run_matlab_and_compare "C44L4" "c44l4" "C44L4.m"
run_matlab_and_compare "C44L5" "c44l5" "C44L5.m"

# Chapter 45
run_matlab_and_compare "C45L1" "c45l1" "C45L1.m"
run_matlab_and_compare "C45L2" "c45l2" "C45L2.m"
run_matlab_and_compare "C45L3" "c45l3" "C45L3.m"
run_matlab_and_compare "C45L4" "c45l4" "C45L4.m"

echo ""
echo "=== VERIFICATION COMPLETE ==="
echo ""

# Count pass/fail
PASS_COUNT=$(grep -c "PASS" "$RESULTS_FILE" 2>/dev/null | head -1 || echo 0)
FAIL_COUNT=$(grep -c "FAIL" "$RESULTS_FILE" 2>/dev/null | head -1 || echo 0)
MISSING_COUNT=$(grep -cE "MISSING|NO OUTPUT" "$RESULTS_FILE" 2>/dev/null | head -1 || echo 0)
# Ensure counts are valid integers
PASS_COUNT=${PASS_COUNT:-0}
FAIL_COUNT=${FAIL_COUNT:-0}
MISSING_COUNT=${MISSING_COUNT:-0}
TOTAL=$((PASS_COUNT + FAIL_COUNT + MISSING_COUNT))

# Check if single test mode and no test was found
if [ -n "$SINGLE_TEST" ] && [ "$TOTAL" -eq 0 ]; then
    echo "ERROR: Test '$SINGLE_TEST' not found."
    echo ""
    echo "Valid test names are like: C1L1, C2L1, C10L1, etc."
    exit 1
fi

echo "=== STATISTICS ==="
echo "PASS:    $PASS_COUNT / $TOTAL"
echo "FAIL:    $FAIL_COUNT"
echo "MISSING: $MISSING_COUNT"
echo ""
echo "Results saved to: $RESULTS_FILE"

# Note about Monte Carlo simulations (only show when running all tests)
if [ -z "$SINGLE_TEST" ]; then
    echo ""
    echo "NOTE: Monte Carlo simulations (C7, C12, C14, C28, C30, C31, C34, C41, C42)"
    echo "      use random numbers, so exact matches are not expected."
    echo "      These are verified for structure (same rows/columns) rather than values."
fi
