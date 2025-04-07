use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
pub struct Cli {
    // ... existing CLI code ...
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    // ... existing commands ...
    /// Process a FASTQ file to classify its contents
    ProcessFastq {
        /// Path to the FASTQ file
        #[arg(short, long)]
        fastq: PathBuf,

        /// Sample ID
        #[arg(short, long)]
        sample_id: String,

        /// Path to the output directory
        #[arg(short, long, default_value = "results")]
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
        #[arg(short, long)]
        dir: PathBuf,

        /// Path to the output directory
        #[arg(short, long, default_value = "results")]
        output: PathBuf,
    },
}

/// Main entry point for CLI
pub fn run_cli(cli: Cli) -> Result<(), Box<dyn std::error::Error>> {
    // ... existing CLI handling ...

    match cli.command {
        // ... existing command handlers ...
        Commands::ProcessFastq {
            fastq,
            sample_id,
            output,
            min_quality,
            min_length,
        } => {
            // Create QC parameters
            let qc_params = QualityControlParams {
                min_avg_quality: min_quality,
                min_length,
                trim_quality: 15,   // Default
                max_n_percent: 5.0, // Default
            };

            // Create FASTQ processor
            let mut processor = FastqProcessor::new(
                cli.db_path,
                cli.cache_dir,
                cli.threads,
                31,   // Default macro_k
                21,   // Default meso_k
                1000, // Default sketch_size
                Some(qc_params),
                cli.api_key,
            )?;

            // Initialize classifier
            processor.init_classifier()?;

            // Process FASTQ file
            let results = processor.process_file(fastq, &sample_id, output)?;

            // Generate and print report
            let report = generate_report(&results)?;
            println!("{}", report);
        }

        Commands::ProcessDir { dir, output } => {
            // Create FASTQ processor
            let mut processor = FastqProcessor::new(
                cli.db_path,
                cli.cache_dir,
                cli.threads,
                31,   // Default macro_k
                21,   // Default meso_k
                1000, // Default sketch_size
                None, // Default QC parameters
                cli.api_key,
            )?;

            // Initialize classifier
            processor.init_classifier()?;

            // Find all FASTQ files in the directory
            let mut fastq_files = Vec::new();
            for entry in std::fs::read_dir(dir)? {
                let entry = entry?;
                let path = entry.path();

                if path.is_file()
                    && path.extension().map_or(false, |ext| {
                        ext == "fastq" || ext == "fq" || ext == "fastq.gz" || ext == "fq.gz"
                    })
                {
                    fastq_files.push(path);
                }
            }

            println!("Found {} FASTQ files to process", fastq_files.len());

            // Process each FASTQ file
            for (i, path) in fastq_files.iter().enumerate() {
                let sample_id = path
                    .file_stem()
                    .unwrap_or_default()
                    .to_string_lossy()
                    .replace(".fastq", "")
                    .replace(".fq", "");

                println!(
                    "Processing file {}/{}: {} (Sample ID: {})",
                    i + 1,
                    fastq_files.len(),
                    path.display(),
                    sample_id
                );

                match processor.process_file(path, &sample_id, &output) {
                    Ok(results) => {
                        println!(
                            "Processed successfully. Wrote results to: {}",
                            results.results_file.unwrap_or_default().display()
                        );
                    }
                    Err(e) => {
                        eprintln!("Error processing {}: {}", path.display(), e);
                    }
                }
            }

            println!("Processed {} FASTQ files", fastq_files.len());
        }
    }

    Ok(())
}
