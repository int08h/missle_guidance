//! Missile Guidance Simulations Runner
//!
//! Command-line interface for running simulations from
//! "Tactical and Strategic Missile Guidance"

// Allow approximate constants (1.5708, 3.1416, 6.28, etc.) to match MATLAB's exact values.
// This is intentional for numerical compatibility - changing these to std::f64::consts
// would cause output differences that violate the MATLAB matching requirement.
#![allow(clippy::approx_constant)]

use std::env;
use std::fs;
use std::io::Write;

mod utils;
mod chapters;
mod plotting;

use chapters::*;

/// Save simulation results to a file in ASCII format (matching MATLAB datfil.txt)
pub fn save_data(filename: &str, data: &[Vec<f64>]) -> std::io::Result<()> {
    let file = std::fs::File::create(filename)?;
    let mut writer = std::io::BufWriter::new(file);

    if data.is_empty() {
        return Ok(());
    }

    let n_rows = data[0].len();

    for i in 0..n_rows {
        for (j, col) in data.iter().enumerate() {
            if j > 0 {
                write!(writer, " ")?;
            }
            write!(writer, "{:14.6e}", col[i])?;
        }
        writeln!(writer)?;
    }

    Ok(())
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        print_help();
        return;
    }

    let command = &args[1];

    match command.as_str() {
        "list" => list_simulations(),
        "run" => {
            if args.len() < 3 {
                println!("Usage: missile_guidance run <simulation>");
                println!("Use 'missile_guidance list' to see available simulations");
                return;
            }
            let sim = &args[2];
            let output_dir = if args.len() > 3 { &args[3] } else { "output" };
            run_simulation(sim, output_dir);
        }
        "run-all" => {
            let output_dir = if args.len() > 2 { &args[2] } else { "output" };
            run_all_simulations(output_dir);
        }
        "help" | "--help" | "-h" => print_help(),
        _ => {
            println!("Unknown command: {}", command);
            print_help();
        }
    }
}

fn print_help() {
    println!("Missile Guidance Simulations");
    println!("============================");
    println!("Rust implementation of 'Tactical and Strategic Missile Guidance' simulations\n");
    println!("Usage:");
    println!("  missile_guidance list              - List all available simulations");
    println!("  missile_guidance run <sim> [dir]   - Run a specific simulation");
    println!("  missile_guidance run-all [dir]     - Run all simulations");
    println!("  missile_guidance help              - Show this help message\n");
    println!("Examples:");
    println!("  missile_guidance run c1l1          - Run Chapter 1, Lesson 1");
    println!("  missile_guidance run c2l1 ./out    - Run C2L1, output to ./out");
    println!("  missile_guidance run-all ./results - Run all simulations");
}

fn list_simulations() {
    println!("Available Simulations:");
    println!("=====================\n");

    for (name, desc) in chapters::list_simulations() {
        println!("  {:8} - {}", name, desc);
    }

    println!("\nUse 'missile_guidance run <name>' to run a simulation");
}

fn run_simulation(sim: &str, output_dir: &str) {
    // Create output directory if it doesn't exist
    if let Err(e) = fs::create_dir_all(output_dir) {
        println!("Error creating output directory: {}", e);
        return;
    }

    println!("Running simulation: {}", sim);
    println!("Output directory: {}\n", output_dir);

    match sim.to_lowercase().as_str() {
        // Chapter 1
        "c1l1" => { c1::l1::run_and_save(output_dir).ok(); }
        "c1l2" => { c1::l2::run_and_save(output_dir).ok(); }
        "c1l3" => { c1::l3::run_and_save(output_dir).ok(); }

        // Chapter 2
        "c2l1" => { c2::l1::run_and_save(output_dir).ok(); }
        "c2l2" => { c2::l2::run_and_save(output_dir).ok(); }

        // Chapter 3
        "c3l1" => { c3::l1::run_and_save(output_dir).ok(); }

        // Chapter 4
        "c4l1" => { c4::l1::run_and_save(output_dir).ok(); }
        "c4l2" => { c4::l2::run_and_save(output_dir).ok(); }
        "c4l3" => { c4::l3::run_and_save(output_dir).ok(); }
        "c4l4" => { c4::l4::run_and_save(output_dir).ok(); }
        "c4l5" => { c4::l5::run_and_save(output_dir).ok(); }
        "c4l6" => { c4::l6::run_and_save(output_dir).ok(); }
        "c4l7" => { c4::l7::run_and_save(output_dir).ok(); }
        "c4l8" => { c4::l8::run_and_save(output_dir).ok(); }

        // Chapter 5
        "c5l1" => { c5::l1::run_and_save(output_dir).ok(); }
        "c5l2" => { c5::l2::run_and_save(output_dir).ok(); }
        "c5l3" => { c5::l3::run_and_save(output_dir).ok(); }

        // Chapter 6
        "c6l1" => { c6::l1::run_and_save(output_dir).ok(); }
        "c6l2" => { c6::l2::run_and_save(output_dir).ok(); }
        "c6l3" => { c6::l3::run_and_save(output_dir).ok(); }
        "c6l4" => { c6::l4::run_and_save(output_dir).ok(); }

        // Chapter 7
        "c7l1" => { c7::l1::run_and_save(output_dir).ok(); }
        "c7l2" => { c7::l2::run_and_save(output_dir).ok(); }
        "c7l3" => { c7::l3::run_and_save(output_dir).ok(); }
        "c7l4" => { c7::l4::run_and_save(output_dir).ok(); }

        // Chapter 8
        "c8l1" => { c8::l1::run_and_save(output_dir).ok(); }
        "c8l2" => { c8::l2::run_and_save(output_dir).ok(); }
        "c8l3" => { c8::l3::run_and_save(output_dir).ok(); }

        // Chapter 9
        "c9l1" => { c9::l1::run_and_save(output_dir).ok(); }
        "c9l2" => { c9::l2::run_and_save(output_dir).ok(); }
        "c9l3" => { c9::l3::run_and_save(output_dir).ok(); }
        "c9l4" => { c9::l4::run_and_save(output_dir).ok(); }
        "c9l5" => { c9::l5::run_and_save(output_dir).ok(); }

        // Chapter 10
        "c10l1" => { c10::l1::run_and_save(output_dir).ok(); }

        // Chapter 11
        "c11l1" => { c11::l1::run_and_save(output_dir).ok(); }
        "c11l2" => { c11::l2::run_and_save(output_dir).ok(); }

        // Chapter 13
        "c13l1" => { c13::l1::run_and_save(output_dir).ok(); }

        // Chapter 12
        "c12l1" => { c12::l1::run_and_save(output_dir).ok(); }
        "c12l2" => { c12::l2::run_and_save(output_dir).ok(); }
        "c12l3" => { c12::l3::run_and_save(output_dir).ok(); }

        // Chapter 14
        "c14l1" => { c14::l1::run_and_save(output_dir).ok(); }
        "c14l2" => { c14::l2::run_and_save(output_dir).ok(); }

        // Chapter 15
        "c15l1" => { c15::l1::run_and_save(output_dir).ok(); }
        "c15l2" => { c15::l2::run_and_save(output_dir).ok(); }
        "c15l3" => { c15::l3::run_and_save(output_dir).ok(); }
        "c15l4" => { c15::l4::run_and_save(output_dir).ok(); }
        "c15l5" => { c15::l5::run_and_save(output_dir).ok(); }
        "c15l6" => { c15::l6::run_and_save(output_dir).ok(); }
        "c15l7" => { c15::l7::run_and_save(output_dir).ok(); }
        "c15l8" => { c15::l8::run_and_save(output_dir).ok(); }

        // Chapter 16
        "c16l1" => { c16::l1::run_and_save(output_dir).ok(); }
        "c16l2" => { c16::l2::run_and_save(output_dir).ok(); }

        // Chapter 17
        "c17l1" => { c17::l1::run_and_save(output_dir).ok(); }
        "c17l2" => { c17::l2::run_and_save(output_dir).ok(); }
        "c17l3" => { c17::l3::run_and_save(output_dir).ok(); }
        "c17l4" => { c17::l4::run_and_save(output_dir).ok(); }
        "c17l5" => { c17::l5::run_and_save(output_dir).ok(); }
        "c17l6" => { c17::l6::run_and_save(output_dir).ok(); }

        // Chapter 18
        "c18l1" => { c18::l1::run_and_save(output_dir).ok(); }
        "c18l2" => { c18::l2::run_and_save(output_dir).ok(); }

        // Chapter 19
        "c19l1" => { c19::l1::run_and_save(output_dir).ok(); }
        "c19l2" => { c19::l2::run_and_save(output_dir).ok(); }

        // Chapter 20
        "c20l1" => { c20::l1::run_and_save(output_dir).ok(); }
        "c20l2" => { c20::l2::run_and_save(output_dir).ok(); }
        "c20l3" => { c20::l3::run_and_save(output_dir).ok(); }
        "c20l4" => { c20::l4::run_and_save(output_dir).ok(); }
        "c20l5" => { c20::l5::run_and_save(output_dir).ok(); }
        "c20l6" => { c20::l6::run_and_save(output_dir).ok(); }
        "c20l7" => { c20::l7::run_and_save(output_dir).ok(); }
        "c20l8" => { c20::l8::run_and_save(output_dir).ok(); }
        "c20l9" => { c20::l9::run_and_save(output_dir).ok(); }

        // Chapter 21
        "c21l1" => { c21::l1::run_and_save(output_dir).ok(); }
        "c21l2" => { c21::l2::run_and_save(output_dir).ok(); }

        // Chapter 22
        "c22l1" => { c22::l1::run_and_save(output_dir).ok(); }
        "c22l2" => { c22::l2::run_and_save(output_dir).ok(); }
        "c22l3" => { c22::l3::run_and_save(output_dir).ok(); }
        "c22l4" => { c22::l4::run_and_save(output_dir).ok(); }
        "c22l5" => { c22::l5::run_and_save(output_dir).ok(); }

        // Chapter 23
        "c23l1" => { c23::l1::run_and_save(output_dir).ok(); }
        "c23l2" => { c23::l2::run_and_save(output_dir).ok(); }
        "c23l3" => { c23::l3::run_and_save(output_dir).ok(); }
        "c23l4" => { c23::l4::run_and_save(output_dir).ok(); }

        // Chapter 24
        "c24l1" => { c24::l1::run_and_save(output_dir).ok(); }
        "c24l2" => { c24::l2::run_and_save(output_dir).ok(); }

        // Chapter 25
        "c25l1" => { c25::l1::run_and_save(output_dir).ok(); }
        "c25l2" => { c25::l2::run_and_save(output_dir).ok(); }
        "c25l3" => { c25::l3::run_and_save(output_dir).ok(); }

        // Chapter 26
        "c26l1" => { c26::l1::run_and_save(output_dir).ok(); }
        "c26l2" => { c26::l2::run_and_save(output_dir).ok(); }
        "c26l3" => { c26::l3::run_and_save(output_dir).ok(); }
        "c26l4" => { c26::l4::run_and_save(output_dir).ok(); }
        "c26l5" => { c26::l5::run_and_save(output_dir).ok(); }
        "c26l6" => { c26::l6::run_and_save(output_dir).ok(); }
        "c26l7" => { c26::l7::run_and_save(output_dir).ok(); }
        "c26l8" => { c26::l8::run_and_save(output_dir).ok(); }
        "c26l9" => { c26::l9::run_and_save(output_dir).ok(); }
        "c26l10" => { c26::l10::run_and_save(output_dir).ok(); }

        // Chapter 27
        "c27l1" => { c27::l1::run_and_save(output_dir).ok(); }
        "c27l2" => { c27::l2::run_and_save(output_dir).ok(); }
        "c27l3" => { c27::l3::run_and_save(output_dir).ok(); }
        "c27l4" => { c27::l4::run_and_save(output_dir).ok(); }
        "c27l5" => { c27::l5::run_and_save(output_dir).ok(); }
        "c27l6" => { c27::l6::run_and_save(output_dir).ok(); }
        "c27l7" => { c27::l7::run_and_save(output_dir).ok(); }

        // Chapter 28
        "c28l1" => { c28::l1::run_and_save(output_dir).ok(); }
        "c28l2" => { c28::l2::run_and_save(output_dir).ok(); }
        "c28l3" => { c28::l3::run_and_save(output_dir).ok(); }

        // Chapter 29
        "c29l1" => { c29::l1::run_and_save(output_dir).ok(); }
        "c29l2" => { c29::l2::run_and_save(output_dir).ok(); }
        "c29l3" => { c29::l3::run_and_save(output_dir).ok(); }
        "c29l4" => { c29::l4::run_and_save(output_dir).ok(); }
        "c29l5" => { c29::l5::run_and_save(output_dir).ok(); }

        // Chapter 30
        "c30l1" => { c30::l1::run_and_save(output_dir).ok(); }
        "c30l2" => { c30::l2::run_and_save(output_dir).ok(); }
        "c30l3" => { c30::l3::run_and_save(output_dir).ok(); }

        // Chapter 31
        "c31l1" => { c31::l1::run_and_save(output_dir).ok(); }

        // Chapter 32
        "c32l1" => { c32::l1::run_and_save(output_dir).ok(); }
        "c32l2" => { c32::l2::run_and_save(output_dir).ok(); }
        "c32l3" => { c32::l3::run_and_save(output_dir).ok(); }
        "c32l4" => { c32::l4::run_and_save(output_dir).ok(); }

        // Chapter 33
        "c33l1" => { c33::l1::run_and_save(output_dir).ok(); }

        // Chapter 34
        "c34l1" => { c34::l1::run_and_save(output_dir).ok(); }
        "c34l2" => { c34::l2::run_and_save(output_dir).ok(); }
        "c34l3" => { c34::l3::run_and_save(output_dir).ok(); }

        // Chapter 35
        "c35l1" => { c35::l1::run_and_save(output_dir).ok(); }
        "c35l2" => { c35::l2::run_and_save(output_dir).ok(); }
        "c35l3" => { c35::l3::run_and_save(output_dir).ok(); }
        "c35l4" => { c35::l4::run_and_save(output_dir).ok(); }
        "c35l5" => { c35::l5::run_and_save(output_dir).ok(); }
        "c35l6" => { c35::l6::run_and_save(output_dir).ok(); }

        // Chapter 36
        "c36l1" => { c36::l1::run_and_save(output_dir).ok(); }
        "c36l2" => { c36::l2::run_and_save(output_dir).ok(); }

        // Chapter 37
        "c37l1" => { c37::l1::run_and_save(output_dir).ok(); }
        "c37l2" => { c37::l2::run_and_save(output_dir).ok(); }
        "c37l3" => { c37::l3::run_and_save(output_dir).ok(); }

        // Chapter 38
        "c38l1" => { c38::l1::run_and_save(output_dir).ok(); }
        "c38l2" => { c38::l2::run_and_save(output_dir).ok(); }
        "c38l3" => { c38::l3::run_and_save(output_dir).ok(); }
        "c38l4" => { c38::l4::run_and_save(output_dir).ok(); }

        // Chapter 39
        "c39l1" => { c39::l1::run_and_save(output_dir).ok(); }
        "c39l2" => { c39::l2::run_and_save(output_dir).ok(); }
        "c39l3" => { c39::l3::run_and_save(output_dir).ok(); }

        // Chapter 40
        "c40l1" => { c40::l1::run_and_save(output_dir).ok(); }
        "c40l2" => { c40::l2::run_and_save(output_dir).ok(); }
        "c40l3" => { c40::l3::run_and_save(output_dir).ok(); }
        "c40l4" => { c40::l4::run_and_save(output_dir).ok(); }
        "c40l5" => { c40::l5::run_and_save(output_dir).ok(); }

        // Chapter 41
        "c41l1" => { c41::l1::run_and_save(output_dir).ok(); }
        "c41l2" => { c41::l2::run_and_save(output_dir).ok(); }

        // Chapter 42
        "c42l1" => { c42::l1::run_and_save(output_dir).ok(); }

        // Chapter 43
        "c43l1" => { c43::l1::run_and_save(output_dir).ok(); }
        "c43l2" => { c43::l2::run_and_save(output_dir).ok(); }

        // Chapter 44
        "c44l1" => { c44::l1::run_and_save(output_dir).ok(); }
        "c44l2" => { c44::l2::run_and_save(output_dir).ok(); }
        "c44l3" => { c44::l3::run_and_save(output_dir).ok(); }
        "c44l4" => { c44::l4::run_and_save(output_dir).ok(); }
        "c44l5" => { c44::l5::run_and_save(output_dir).ok(); }

        // Chapter 45
        "c45l1" => { c45::l1::run_and_save(output_dir).ok(); }
        "c45l2" => { c45::l2::run_and_save(output_dir).ok(); }
        "c45l3" => { c45::l3::run_and_save(output_dir).ok(); }
        "c45l4" => { c45::l4::run_and_save(output_dir).ok(); }

        _ => {
            println!("Unknown simulation: {}", sim);
            println!("Use 'missile_guidance list' to see available simulations");
        }
    }
}

fn run_all_simulations(output_dir: &str) {
    // Create output directory if it doesn't exist
    if let Err(e) = fs::create_dir_all(output_dir) {
        println!("Error creating output directory: {}", e);
        return;
    }

    println!("Running all simulations...");
    println!("Output directory: {}\n", output_dir);

    let sims = [
        "c1l1", "c1l2", "c1l3",
        "c2l1", "c2l2",
        "c3l1",
        "c4l1", "c4l2", "c4l3", "c4l4", "c4l5", "c4l6", "c4l7", "c4l8",
        "c5l1", "c5l2", "c5l3",
        "c6l1", "c6l2", "c6l3", "c6l4",
        "c7l1", "c7l2", "c7l3", "c7l4",
        "c8l1", "c8l2", "c8l3",
        "c9l1", "c9l2", "c9l3", "c9l4", "c9l5",
        "c10l1",
        "c11l1", "c11l2",
        "c12l1", "c12l2", "c12l3",
        "c13l1",
        "c14l1", "c14l2",
        "c15l1", "c15l2", "c15l3", "c15l4", "c15l5", "c15l6", "c15l7", "c15l8",
        "c16l1", "c16l2",
        "c17l1", "c17l2", "c17l3", "c17l4", "c17l5", "c17l6",
        "c18l1", "c18l2",
        "c19l1", "c19l2",
        "c20l1", "c20l2", "c20l3", "c20l4", "c20l5", "c20l6", "c20l7", "c20l8", "c20l9",
        "c21l1", "c21l2",
        "c22l1", "c22l2", "c22l3", "c22l4", "c22l5",
        "c23l1", "c23l2", "c23l3", "c23l4",
        "c24l1", "c24l2",
        "c25l1", "c25l2", "c25l3",
        "c26l1", "c26l2", "c26l3", "c26l4", "c26l5", "c26l6", "c26l7", "c26l8", "c26l9", "c26l10",
        "c27l1", "c27l2", "c27l3", "c27l4", "c27l5", "c27l6", "c27l7",
        "c28l1", "c28l2", "c28l3",
        "c29l1", "c29l2", "c29l3", "c29l4", "c29l5",
        "c30l1", "c30l2", "c30l3",
        "c31l1",
        "c32l1", "c32l2", "c32l3", "c32l4",
        "c33l1",
        "c34l1", "c34l2", "c34l3",
        "c35l1", "c35l2", "c35l3", "c35l4", "c35l5", "c35l6",
        "c36l1", "c36l2",
        "c37l1", "c37l2", "c37l3",
        "c38l1", "c38l2", "c38l3", "c38l4",
        "c39l1", "c39l2", "c39l3",
        "c40l1", "c40l2", "c40l3", "c40l4", "c40l5",
        "c41l1", "c41l2",
        "c42l1",
        "c43l1", "c43l2",
        "c44l1", "c44l2", "c44l3", "c44l4", "c44l5",
        "c45l1", "c45l2", "c45l3", "c45l4",
    ];

    for sim in &sims {
        println!("\n--- Running {} ---", sim);
        run_simulation(sim, output_dir);
    }

    println!("\n=== All simulations complete ===");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c1l1_simulation() {
        let results = c1::l1::run();
        assert!(!results.time.is_empty());
    }

    #[test]
    fn test_c2l1_simulation() {
        let results = c2::l1::run();
        assert!(!results.time.is_empty());
    }
}
