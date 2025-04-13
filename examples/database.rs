use std::error::Error;
use std::path::Path;
use strain_ahsp::database::DatabaseManager;

fn main() -> Result<(), Box<dyn Error>> {
    // Initialize database manager
    let mut manager = DatabaseManager::new(
        "ahsp_db",
        "genome_cache",
        31,   // macro_k
        21,   // meso_k
        1000, // sketch_size
        8,    // threads
        None, // api_key
    )?;

    // Check if database is empty
    if manager.is_empty()? {
        println!("Initializing database with reference genomes...");

        // Download and process reference genomes
        let added = manager.search_and_add_references("Escherichia coli", 5)?;
        println!("Added {} references", added.len());
    }

    // Get all signatures for classification
    let references = manager.database.get_all_signatures()?;

    // Create adaptive classifier
    let classifier = AdaptiveClassifier::new(references, None, None)
        .map_err(|e| Box::new(e) as Box<dyn Error>)?;

    // Now you can classify metagenomic samples
    println!("Ready to classify samples");

    Ok(())
}
