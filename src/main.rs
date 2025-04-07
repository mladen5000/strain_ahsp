//! Main entry point for the metagenomics_deseq2 application.
//!
//! This application likely performs differential abundance analysis
//! on metagenomic data, potentially using methods inspired by DESeq2.
//! It might involve steps like:
//! 1. Reading sequencing data (FASTQ).
//! 2. Generating k-mer counts or using sketching techniques.
//! 3. Building a count table (features vs. samples).
//! 4. Normalizing counts.
//! 5. Performing statistical tests for differential abundance.
//! 6. Outputting results.

// Modules defined within the project
mod adaptive;
mod bio;
mod cli;
mod config;
mod count_table;
mod database;
mod io;
mod midas_db;
mod normalization;
mod pipeline;
mod sketch;
mod stats;
mod strain_method;
mod utils;
use cli::{run_cli, Cli};

// External Crate Imports
use anyhow::Result;
use clap::Parser;
use log::{error, info}; // Using log crate for logging // Using anyhow for easier error handling

// Local Imports
// Example: use crate::io::fastq::read_fastq_files;
// Example: use crate::count_table::CountTable;
// Example: use crate::normalization::normalize_counts;
// Example: use crate::stats::perform_differential_analysis;

/// Define command-line arguments using clap.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input FASTQ file(s) or directory containing FASTQ files.
    #[arg(short, long, required = true)]
    input: Vec<String>,

    /// Output file path for the results (e.g., CSV).
    #[arg(short, long, required = true)]
    output: String,

    /// Metadata file describing samples and conditions (e.g., CSV).
    #[arg(short, long)]
    metadata: Option<String>,

    /// K-mer size for analysis.
    #[arg(short = 'k', long, default_value_t = 31)]
    kmer_size: usize,

    /// Number of threads to use.
    #[arg(short = 't', long, default_value_t = 1)]
    threads: usize,

    /// Normalization method to apply (e.g., "median-of-ratios", "tpm", "none").
    #[arg(long, default_value = "median-of-ratios")]
    normalization: String,
    // Add other relevant arguments:
    // - Sketch size if using MinHash/sketching
    // - Minimum count threshold
    // - Significance level (alpha)
    // - MIDAS DB path if used
    // - Strain analysis parameters
    // ...
}

/// Main function: parses arguments and orchestrates the analysis workflow.
fn main() -> Result<()> {
    // Initialize logging (e.g., using env_logger)
    env_logger::init();

    // Parse command-line arguments
    let args = Args::parse();
    info!("Starting analysis with arguments: {:?}", args);

    // Configure Rayon thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()?;
    info!("Using {} threads.", args.threads);

    // --- Core Workflow ---
    // The following steps are placeholders and need specific implementation.

    // 1. Read Input Data (e.g., FASTQ files)
    info!("Reading input files: {:?}", args.input);
    // let sequences = io::fastq::read_sequences(&args.input)?;
    // TODO: Implement sequence reading, potentially in parallel.

    // 2. Process Sequences (e.g., K-mer counting or Sketching)
    info!("Processing sequences (k={})...", args.kmer_size);
    // let counts_or_sketches = bio::kmers::process_sequences(sequences, args.kmer_size)?;
    // or
    // let sketches = sketch::minhash::calculate_sketches(sequences, args.kmer_size, sketch_size)?;
    // TODO: Implement k-mer counting or sketching logic.

    // 3. Build Count Table
    info!("Building count table...");
    // let count_table = count_table::build_table(counts_or_sketches)?;
    // TODO: Implement logic to aggregate counts/sketches into a table (features x samples).
    // TODO: Potentially load metadata here if needed for table structure.

    // 4. Normalize Count Table
    info!(
        "Normalizing counts using '{}' method...",
        args.normalization
    );
    // let normalized_table = normalization::normalize(&count_table, &args.normalization)?;
    // TODO: Implement different normalization strategies (Median-of-Ratios like DESeq2, TPM, etc.).

    // 5. Perform Statistical Analysis (Differential Abundance)
    info!("Performing differential abundance analysis...");
    // let results = stats::run_deseq2_like_analysis(&normalized_table, &args.metadata)?;
    // or potentially using Bayesian methods:
    // let results = stats::bayesian::run_analysis(&normalized_table, &args.metadata)?;
    // TODO: Implement the core statistical testing logic. This is complex.
    // TODO: May need to parse metadata file here.

    // 6. Strain-Level Analysis (Optional)
    // if strain_analysis_requested {
    //     info!("Performing strain-level analysis...");
    //     let strain_results = strain_method::analyze_strains(/* relevant data */)?;
    //     // TODO: Implement strain analysis logic.
    // }

    // 7. Write Output
    info!("Writing results to {}...", args.output);
    // io::write_results(&results, &args.output)?;
    // TODO: Implement result writing (e.g., to CSV).

    info!("Analysis finished successfully.");
    // Initialize logger
    env_logger::init();

    // Parse command line arguments
    let cli = Cli::parse();

    // Run CLI
    run_cli(cli)?;

    Ok(())
}

// Unit tests module
#[cfg(test)]
mod tests {
    // TODO: Add unit tests for functions in main.rs if any logic is implemented here.
    // Example:
    // #[test]
    // fn test_example() {
    //     assert_eq!(2 + 2, 4);
    // }
}
