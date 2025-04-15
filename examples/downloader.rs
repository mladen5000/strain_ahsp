// examples/populate_database.rs

use clap::Parser;
use std::path::PathBuf;
use std::process::exit;
use std::time::Instant;

// Assuming your library crate is named 'qc_sketch_db' (adjust if different)
// You might need to adjust the path based on your project structure if database.rs isn't in src/lib.rs
use strain_ahsp::database::DatabaseManager; // Use your actual crate name

/// Command-line arguments for the database population tool
#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// NCBI search query (e.g., "Escherichia coli", "Bacillus subtilis")
    #[clap(short, long, value_parser)]
    query: String,

    /// Directory to store the signature database
    #[clap(short, long, value_parser, default_value = "signature_db")]
    db_dir: PathBuf,

    /// Directory to cache downloaded genome files
    #[clap(short, long, value_parser, default_value = "ncbi_cache")]
    cache_dir: PathBuf,

    /// Maximum number of genomes to download and process
    #[clap(short, long, value_parser, default_value_t = 10)]
    max_results: usize,

    /// Primary k-mer size for the SignatureBuilder
    #[clap(long, value_parser, default_value_t = 31)]
    kmer_size: usize,

    /// Sketch size (or scaling factor) for the SignatureBuilder
    #[clap(long, value_parser, default_value_t = 1000)]
    sketch_size: usize,

    /// Optional NCBI API key
    #[clap(long, value_parser)]
    api_key: Option<String>,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Initialize logging (shows info!, warn!, error! messages)
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    // Parse command-line arguments
    let args = Args::parse();

    log::info!("Starting database population process...");
    log::info!("Search Query: {}", args.query);
    log::info!("Database Directory: {}", args.db_dir.display());
    log::info!("Cache Directory: {}", args.cache_dir.display());
    log::info!("Max Results: {}", args.max_results);
    log::info!("Signature Builder k: {}", args.kmer_size);
    log::info!("Signature Builder size: {}", args.sketch_size);
    log::info!("NCBI API Key Provided: {}", args.api_key.is_some());

    // --- Create directories if they don't exist ---
    if let Err(e) = std::fs::create_dir_all(&args.db_dir) {
        log::error!(
            "Failed to create database directory '{}': {}",
            args.db_dir.display(),
            e
        );
        exit(1);
    }
    if let Err(e) = std::fs::create_dir_all(&args.cache_dir) {
        log::error!(
            "Failed to create cache directory '{}': {}",
            args.cache_dir.display(),
            e
        );
        exit(1);
    }

    // --- Initialize the Database Manager ---
    // NOTE: The DatabaseManager::new implementation provided in the prompt
    //       hardcodes SignatureBuilder::new(31, 21, 1000, 1).
    //       However, the *signature* takes builder_kmer_size and builder_sketch_size.
    //       This example follows the *signature*, passing the command-line args.
    //       Ensure your actual DatabaseManager::new implementation uses these parameters
    //       or adjust this example accordingly.
    log::info!("Initializing Database Manager...");
    let mut manager = match DatabaseManager::new(
        &args.db_dir,
        &args.cache_dir,
        args.kmer_size,   // Pass k-mer size from args
        args.sketch_size, // Pass sketch size from args
        args.api_key.clone(),
    ) {
        Ok(m) => m,
        Err(e) => {
            log::error!("Failed to initialize Database Manager: {}", e);
            exit(1);
        }
    };

    // --- Search, Download, and Process References ---
    log::info!(
        "Searching NCBI for '{}' and processing up to {} results...",
        args.query,
        args.max_results
    );
    let start_time = Instant::now();

    match manager.search_and_add_references(&args.query, args.max_results) {
        Ok(added_ids) => {
            let duration = start_time.elapsed();
            if added_ids.is_empty() {
                log::info!(
                    "Completed successfully, but no new genomes matching the criteria were found or processed. (Duration: {:.2?})",
                    duration
                );
            } else {
                log::info!(
                    "Successfully added {} new signatures to the database. (Duration: {:.2?})",
                    added_ids.len(),
                    duration
                );
                log::debug!("Added accessions: {:?}", added_ids);
            }
        }
        Err(e) => {
            let duration = start_time.elapsed();
            log::error!(
                "An error occurred during the process: {} (Duration: {:.2?})",
                e,
                duration
            );
            // Depending on the error, you might want different exit codes
            exit(1);
        }
    }

    log::info!("Database population process finished.");
    Ok(())
}
