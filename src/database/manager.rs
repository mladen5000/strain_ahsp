use clap::{Parser, Subcommand};
use std::path::PathBuf;

use crate::database::DatabaseManager;
use log::{error, info, warn}; // Added log imports

#[derive(Parser, Debug)] // Added Debug
#[command(author, version, about, long_about = None)]
pub struct Cli {
    /// Path to the signature database directory
    #[arg(long, value_name = "DIR", default_value = "ahsp_db")] // Use long flags for clarity
    pub db_path: PathBuf,

    /// Path to the genome cache directory
    #[arg(long, value_name = "DIR", default_value = "genome_cache")]
    pub cache_dir: PathBuf,

    /// NCBI API key (optional)
    #[arg(long)] // Use long flag
    pub api_key: Option<String>,

    // threads is a global option but likely not needed for DatabaseManager::new itself
    // Rayon will pick it up, or you can configure Rayon globally if needed.
    /// Number of threads for parallel operations (e.g., downloads)
    #[arg(short, long, default_value_t = 4)]
    pub threads: usize,

    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)] // Added Debug
pub enum Commands {
    /// Initialize the database with reference genomes
    Init {
        /// Query to search for reference genomes on NCBI
        #[arg(short, long, required = true)] // Mark as required
        query: String,

        /// Maximum number of reference genomes to download
        #[arg(long, default_value_t = 20)] // Use long flag
        max_refs: usize,

        /// K-mer size for the primary signature level (e.g., macro)
        #[arg(long, default_value_t = 31)]
        kmer_size: usize, // Renamed from macro_k for clarity in `new` call

        // meso_k is specific to the builder's internal logic, not needed for DatabaseManager::new
        // If your actual builder *requires* it in ::new, adjust DatabaseManager::new signature
        // #[arg(long, default_value_t = 21)]
        // meso_k: usize,
        /// Sketch size (number of hashes) for MinHash signatures
        #[arg(long, default_value_t = 1000)]
        sketch_size: usize,
    },

    /// Add new reference genomes to the database
    AddReferences {
        /// Query to search for reference genomes on NCBI
        #[arg(short, long, required = true)] // Mark as required
        query: String,

        /// Maximum number of reference genomes to download
        #[arg(long, default_value_t = 10)] // Use long flag
        max_refs: usize,
        // Note: Uses default kmer/sketch sizes when adding references later
    },

    /// List all reference genomes currently in the database
    ListReferences,

    /// Search for signatures matching a taxonomy term (name or ID)
    Search {
        /// Taxonomy term (e.g., 'Escherichia coli', '562') to search for
        #[arg(short, long, required = true)] // Mark as required
        term: String,
    },
}

/// Main entry point for database management CLI
pub fn run_database_cli(cli: Cli) -> Result<(), Box<dyn std::error::Error>> {
    // Initialize logger (only once)
    // Consider using a more robust logger setup like fern or tracing
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    // Configure Rayon thread pool if explicit control is needed
    // rayon::ThreadPoolBuilder::new().num_threads(cli.threads).build_global().unwrap();
    // info!("Set Rayon global thread pool to {} threads", cli.threads);
    // Otherwise, Rayon typically uses the number of logical cores by default.

    match cli.command {
        Commands::Init {
            query,
            max_refs,
            kmer_size, // Use the renamed argument
            // meso_k is not passed here
            sketch_size,
        } => {
            info!("Initializing database...");
            // Create database manager with parameters from Init command
            // Assuming DatabaseManager::new(db, cache, k, sketch, api_key) -> 5 args
            let mut manager = DatabaseManager::new(
                &cli.db_path,        // Pass as reference
                &cli.cache_dir,      // Pass as reference
                kmer_size,           // Pass k-mer size from command
                sketch_size,         // Pass sketch size from command
                cli.api_key.clone(), // Clone Option<String>
            )?;
            info!(
                "DatabaseManager created with k={}, sketch_size={}",
                kmer_size, sketch_size
            );

            // Check if database is already initialized
            // Note: DatabaseManager::is_empty checks signature count, not just existence of DB files
            if !manager.is_empty()? {
                warn!(
                    "Database at '{}' already contains signatures. Use 'add-references' command to add more.",
                    cli.db_path.display()
                );
                // Optionally, add a --force flag to Init to allow overwriting/clearing
                return Ok(());
            }

            // Initialize database with references
            info!(
                "Populating database with initial references for query: '{}' (max: {})",
                query, max_refs
            );
            let added_ids = manager.search_and_add_references(&query, max_refs)?;

            if added_ids.is_empty() {
                info!("No reference signatures were added (query might have yielded no results or downloads failed).");
            } else {
                info!(
                    "Successfully added {} reference signatures:",
                    added_ids.len()
                );
                for id in added_ids {
                    println!("  - {}", id); // Keep simple output for user
                }
            }
            info!("Database initialization complete.");
        }

        Commands::AddReferences { query, max_refs } => {
            info!("Adding references to existing database...");
            // Create database manager with default signature parameters
            // Assuming DatabaseManager::new(db, cache, k, sketch, api_key) -> 5 args
            let mut manager = DatabaseManager::new(
                &cli.db_path,
                &cli.cache_dir,
                31,   // Default k-mer size for adding later
                1000, // Default sketch size for adding later
                cli.api_key.clone(),
            )?;
            info!("DatabaseManager created with default signature parameters (k=31, sketch=1000)");

            // Add references
            info!(
                "Searching and adding references for query: '{}' (max: {})",
                query, max_refs
            );
            let added_ids = manager.search_and_add_references(&query, max_refs)?;

            if added_ids.is_empty() {
                info!("No new reference signatures were added.");
            } else {
                info!(
                    "Successfully added {} new reference signatures:",
                    added_ids.len()
                );
                for id in added_ids {
                    println!("  - {}", id);
                }
            }
        }

        Commands::ListReferences => {
            info!("Listing references from database...");
            // Create database manager - signature params don't matter for listing
            let manager = DatabaseManager::new(
                &cli.db_path,
                &cli.cache_dir,
                31,   // Default k-mer size (arbitrary for this command)
                1000, // Default sketch size (arbitrary for this command)
                cli.api_key.clone(),
            )?;

            // List all references
            let signatures = manager.database.get_all_signatures()?;

            if signatures.is_empty() {
                println!(
                    "Database at '{}' contains 0 reference signatures.",
                    cli.db_path.display()
                );
            } else {
                println!(
                    "Database at '{}' contains {} reference signatures:",
                    cli.db_path.display(),
                    signatures.len()
                );
                // Sort signatures by ID for consistent output
                let mut sorted_signatures = signatures;
                sorted_signatures.sort_by(|a, b| a.taxon_id.cmp(&b.taxon_id));

                for sig in sorted_signatures {
                    // Attempt to get the species name (often the last element)
                    let species_name = sig
                        .lineage
                        .last()
                        .cloned()
                        .unwrap_or_else(|| "Unknown Species".to_string());
                    println!("  - {} ({})", sig.taxon_id, species_name);
                }
            }
        }

        Commands::Search { term } => {
            info!("Searching database for term: '{}'", term);
            // Create database manager - signature params don't matter for searching
            let manager = DatabaseManager::new(
                &cli.db_path,
                &cli.cache_dir,
                31,   // Default k-mer size (arbitrary for this command)
                1000, // Default sketch size (arbitrary for this command)
                cli.api_key.clone(),
            )?;

            // Search for signatures
            let results = manager.database.search_by_taxonomy(&term)?;

            if results.is_empty() {
                println!(
                    "No signatures found matching term '{}' in database '{}'.",
                    term,
                    cli.db_path.display()
                );
            } else {
                println!("Found {} matches for term '{}':", results.len(), term);
                // Sort results by ID for consistent output
                let mut sorted_results = results;
                sorted_results.sort_by(|a, b| a.taxon_id.cmp(&b.taxon_id));

                for sig in sorted_results {
                    let species_name = sig
                        .lineage
                        .last()
                        .cloned()
                        .unwrap_or_else(|| "Unknown Species".to_string());
                    println!("  - {} ({})", sig.taxon_id, species_name);
                }
            }
        }
    }

    Ok(())
}
