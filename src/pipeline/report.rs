use clap::{Parser, Subcommand};
use log::info;
use std::path::PathBuf;

// Assuming these imports are correct relative to your project structure
use crate::pipeline::{
    processor::{generate_report, QualityControlParams},
    FastqProcessor,
};

#[derive(Parser, Debug)] // Added Debug for easier printing if needed
#[command(author, version, about, long_about = None)]
pub struct Cli {
    /// Path to the signature database directory
    #[arg(long, value_name = "DIR", required = true)] // Made required explicitly
    pub db_path: PathBuf,

    /// Path to the cache directory for downloads
    #[arg(long, value_name = "DIR", required = true)] // Made required explicitly
    pub cache_dir: PathBuf,

    /// Number of threads to use for processing
    #[arg(short, long, default_value_t = 4)] // Set a default value
    pub threads: usize,

    /// NCBI API key (optional)
    #[arg(long)]
    pub api_key: Option<String>,

    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)] // Added Debug
pub enum Commands {
    /// Process a FASTQ file to classify its contents
    ProcessFastq {
        /// Path to the FASTQ file
        #[arg(short, long, value_name = "FILE", required = true)]
        fastq: PathBuf,

        /// Sample ID
        #[arg(short, long, required = true)]
        sample_id: String,

        /// Path to the output directory
        #[arg(short, long, default_value = "results", value_name = "DIR")]
        output: PathBuf,

        /// Minimum average quality score for reads
        #[arg(long, default_value_t = 20.0)]
        min_quality: f64,

        /// Minimum read length after trimming
        #[arg(long, default_value_t = 50)]
        min_length: usize,
    },

    /// Process multiple FASTQ files in a directory
    ProcessDir {
        /// Path to the directory containing FASTQ files
        #[arg(short, long, value_name = "DIR", required = true)]
        dir: PathBuf,

        /// Path to the output directory
        #[arg(short, long, default_value = "results", value_name = "DIR")]
        output: PathBuf,
    },
}

/// Main entry point for CLI
pub fn run_cli(cli: Cli) -> Result<(), Box<dyn std::error::Error>> {
    // Configure logging (example using env_logger) - add if you haven't
    // env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    // Now you can access db_path, cache_dir etc. directly from cli *before* the match

    match cli.command {
        Commands::ProcessFastq {
            fastq,
            sample_id,
            output,
            min_quality,
            min_length,
        } => {
            info!(
                "Processing FASTQ file: {} with Sample ID: {}",
                fastq.display(),
                sample_id
            );

            // Create QC parameters
            let qc_params = QualityControlParams {
                min_avg_quality: min_quality,
                min_length,
                trim_quality: 15,   // Example Default
                max_n_percent: 5.0, // Example Default
            };
            info!("QC Parameters: {:?}", qc_params);

            // Create FASTQ processor using the global args from `cli`
            let mut processor = FastqProcessor::new(
                &cli.db_path,        // Pass reference if needed by constructor
                &cli.cache_dir,      // Pass reference if needed by constructor
                cli.threads,         // Pass value
                31,                  // Default macro_k - consider making these CLI args too?
                21,                  // Default meso_k
                1000,                // Default sketch_size
                Some(qc_params),     // Pass specific QC params for this command
                cli.api_key.clone(), // Clone Option<String> if needed
            )?;
            info!("FastqProcessor created.");

            // Initialize classifier
            processor.init_classifier()?;
            info!("Classifier initialized.");

            // Process FASTQ file
            let results = processor.process_file(&fastq, &sample_id, &output)?; // Pass references
            info!("File processing complete. Results: {:?}", results); // Example log

            // Generate and print report
            // Ensure generate_report takes the correct type from process_file result
            // let report = generate_report(&results)?;
            // println!("{}", report);
            println!("Processing finished. Results summary struct: {:?}", results);
            // Placeholder report
        }

        Commands::ProcessDir { dir, output } => {
            info!(
                "Processing directory: {} into output: {}",
                dir.display(),
                output.display()
            );

            // Create FASTQ processor using the global args from `cli`
            let mut processor = FastqProcessor::new(
                &cli.db_path,        // Pass reference
                &cli.cache_dir,      // Pass reference
                cli.threads,         // Pass value
                31,                  // Default macro_k
                21,                  // Default meso_k
                1000,                // Default sketch_size
                None, // No specific QC parameters for directory processing (uses defaults in processor)
                cli.api_key.clone(), // Clone Option<String>
            )?;
            info!("FastqProcessor created for directory processing.");

            // Initialize classifier
            processor.init_classifier()?;
            info!("Classifier initialized.");

            // Find all FASTQ files in the directory
            let mut fastq_files = Vec::new();
            for entry in std::fs::read_dir(&dir)? {
                // Pass reference to dir
                let entry = entry?;
                let path = entry.path();

                // Improved check for fastq files (case-insensitive extensions)
                if path.is_file() {
                    if let Some(ext) = path.extension().and_then(|e| e.to_str()) {
                        let lower_ext = ext.to_lowercase();
                        if lower_ext == "fastq" || lower_ext == "fq" {
                            fastq_files.push(path);
                        } else if lower_ext == "gz" {
                            // Check the part before .gz
                            if let Some(stem) = path.file_stem() {
                                if let Some(stem_str) = stem.to_str() {
                                    let stem_path = PathBuf::from(stem_str);
                                    if let Some(stem_ext) = stem_path.extension() {
                                        if let Some(ext_str) = stem_ext.to_str() {
                                            let lower_stem_ext = ext_str.to_lowercase();
                                            if lower_stem_ext == "fastq" || lower_stem_ext == "fq" {
                                                fastq_files.push(path);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if fastq_files.is_empty() {
                log::warn!(
                    "No FASTQ files (.fastq, .fq, .fastq.gz, .fq.gz) found in directory: {}",
                    dir.display()
                );
                return Ok(()); // Nothing to do
            }

            println!("Found {} FASTQ files to process.", fastq_files.len());

            // Process each FASTQ file
            for (i, path) in fastq_files.iter().enumerate() {
                // Generate sample ID from file stem more robustly
                let sample_id = path
                    .file_name() // Get full filename first
                    .and_then(|name| name.to_str())
                    .map(|name_str| {
                        // Remove common fastq extensions
                        name_str
                            .trim_end_matches(".gz")
                            .trim_end_matches(".fastq")
                            .trim_end_matches(".fq")
                            .to_string()
                    })
                    .unwrap_or_else(|| format!("sample_{}", i + 1)); // Fallback ID

                println!(
                    "Processing file {}/{}: {} (Sample ID: {})",
                    i + 1,
                    fastq_files.len(),
                    path.display(),
                    sample_id
                );

                // Ensure paths are passed as references
                match processor.process_file(path, &sample_id, &output) {
                    Ok(results) => {
                        println!(
                            "Processed '{}' successfully. Results file: {}",
                            sample_id,
                            // Handle option properly if results_file can be None
                            results
                                .results_file
                                .as_ref()
                                .map(|p| p.display().to_string())
                                .unwrap_or_else(|| "N/A".to_string())
                        );
                    }
                    Err(e) => {
                        eprintln!("Error processing {}: {}", path.display(), e);
                        // Decide if you want to continue processing other files or stop
                        // continue; // Example: continue to next file on error
                    }
                }
            }

            println!("Finished processing {} FASTQ files.", fastq_files.len());
        }
    }

    Ok(())
}
