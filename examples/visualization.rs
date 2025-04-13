use crate::main::FastqProcessor;
fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Process a FASTQ file
    let mut processor = FastqProcessor::new(
        "ahsp_db",
        "genome_cache",
        8,    // threads
        31,   // macro_k
        21,   // meso_k
        1000, // sketch_size
        None, // Default QC parameters
        None, // No API key
    )?;

    // Initialize classifier
    processor.init_classifier()?;

    // Process FASTQ file
    let results = processor.process_file("sample.fastq", "sample1", "results")?;

    // Generate visualizations
    let visualization_files = processor.generate_visualizations(&results, "visualizations")?;

    println!("Generated {} visualizations:", visualization_files.len());
    for file in visualization_files {
        println!("  - {}", file.display());
    }

    Ok(())
}
