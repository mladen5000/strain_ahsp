use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
pub struct Cli {
    /// Path to the signature database
    #[arg(short, long, default_value = "ahsp_db")]
    pub db_path: PathBuf,

    /// Path to the genome cache directory
    #[arg(short, long, default_value = "genome_cache")]
    pub cache_dir: PathBuf,

    /// NCBI API key (optional)
    #[arg(short, long)]
    pub api_key: Option<String>,

    /// Number of threads to use
    #[arg(short, long, default_value_t = 4)]
    pub threads: usize,

    /// Subcommand
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Initialize the database with reference genomes
    Init {
        /// Query to search for reference genomes
        #[arg(short, long)]
        query: String,

        /// Maximum number of references to download
        #[arg(short, long, default_value_t = 20)]
        max_refs: usize,

        /// K-mer size for macro-level signatures
        #[arg(long, default_value_t = 31)]
        macro_k: usize,

        /// K-mer size for meso-level signatures
        #[arg(long, default_value_t = 21)]
        meso_k: usize,

        /// Sketch size for MinHash signatures
        #[arg(long, default_value_t = 1000)]
        sketch_size: usize,
    },

    /// Add new reference genomes to the database
    AddReferences {
        /// Query to search for reference genomes
        #[arg(short, long)]
        query: String,

        /// Maximum number of references to download
        #[arg(short, long, default_value_t = 10)]
        max_refs: usize,
    },

    /// List all reference genomes in the database
    ListReferences,

    /// Search for signatures in the database
    Search {
        /// Taxonomy term to search for
        #[arg(short, long)]
        term: String,
    },
}

/// Main entry point for database management
pub fn run_database_cli(cli: Cli) -> Result<(), Box<dyn std::error::Error>> {
    // Initialize logger
    env_logger::init();

    match cli.command {
        Commands::Init {
            query,
            max_refs,
            macro_k,
            meso_k,
            sketch_size,
        } => {
            // Create database manager
            let mut manager = DatabaseManager::new(
                cli.db_path,
                cli.cache_dir,
                macro_k,
                meso_k,
                sketch_size,
                cli.threads,
                cli.api_key,
            )?;

            // Check if database is already initialized
            if !manager.is_empty()? {
                eprintln!(
                    "Database already contains signatures. Use 'add-references' to add more."
                );
                return Ok(());
            }

            // Initialize database with references
            println!("Initializing database with references for query: {}", query);
            let added = manager.search_and_add_references(&query, max_refs)?;

            println!("Added {} reference signatures:", added.len());
            for id in added {
                println!("  - {}", id);
            }
        }

        Commands::AddReferences { query, max_refs } => {
            // Create database manager
            let mut manager = DatabaseManager::new(
                cli.db_path,
                cli.cache_dir,
                31,   // Default macro_k
                21,   // Default meso_k
                1000, // Default sketch_size
                cli.threads,
                cli.api_key,
            )?;

            // Add references
            println!("Adding references for query: {}", query);
            let added = manager.search_and_add_references(&query, max_refs)?;

            println!("Added {} reference signatures:", added.len());
            for id in added {
                println!("  - {}", id);
            }
        }

        Commands::ListReferences => {
            // Create database manager
            let manager = DatabaseManager::new(
                cli.db_path,
                cli.cache_dir,
                31,   // Default macro_k
                21,   // Default meso_k
                1000, // Default sketch_size
                cli.threads,
                cli.api_key,
            )?;

            // List all references
            let signatures = manager.database.get_all_signatures()?;

            println!(
                "Database contains {} reference signatures:",
                signatures.len()
            );
            for sig in signatures {
                println!(
                    "  - {} ({})",
                    sig.taxon_id,
                    sig.lineage.last().unwrap_or(&"Unknown".to_string())
                );
            }
        }

        Commands::Search { term } => {
            // Create database manager
            let manager = DatabaseManager::new(
                cli.db_path,
                cli.cache_dir,
                31,   // Default macro_k
                21,   // Default meso_k
                1000, // Default sketch_size
                cli.threads,
                cli.api_key,
            )?;

            // Search for signatures
            let results = manager.database.search_by_taxonomy(&term)?;

            println!("Found {} matches for term '{}':", results.len(), term);
            for sig in results {
                println!(
                    "  - {} ({})",
                    sig.taxon_id,
                    sig.lineage.last().unwrap_or(&"Unknown".to_string())
                );
            }
        }
    }

    Ok(())
}
