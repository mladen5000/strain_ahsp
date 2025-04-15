use strain_ahsp::sketch::minhash::{
    AdaptiveClassifier, MultiResolutionSignature, SignatureBuilder,
};
use strain_ahsp::stats::StrainMixtureModel;

use std::error::Error;
use std::path::Path;

fn main() -> Result<(), Box<dyn Error>> {
    // 1. Build reference signatures
    println!("Building reference signatures...");
    let builder = SignatureBuilder::new(31, 21, 1000, 8)?;

    let reference_files = vec![
        (
            Path::new("references/reference1.fasta"),
            "strain1".to_string(),
            vec![
                "Bacteria".to_string(),
                "Proteobacteria".to_string(), /* ... */
            ],
        ),
        (
            Path::new("references/reference2.fasta"),
            "strain2".to_string(),
            vec![
                "Bacteria".to_string(),
                "Firmicutes".to_string(), /* ... */
            ],
        ),
        // More references...
    ];

    let reference_signatures = builder.build_batch(reference_files)?;

    // 2. Create classifier
    println!("Creating classifier...");
    let classifier = AdaptiveClassifier::new(reference_signatures, None, None)?;

    // 3. Process query sample
    println!("Processing query sample...");
    let query_signature = builder.build_from_file(
        "query.fasta",
        "query",
        Vec::new(), // No lineage for query
    )?;

    // 4. Classify at appropriate resolution
    println!("Classifying...");
    let classification = classifier.classify(&query_signature)?;

    println!("Classification result:");
    println!("  Taxon ID: {}", classification.taxon_id);
    println!("  Taxonomic level: {:?}", classification.level);
    println!("  Confidence: {:.2}", classification.confidence);
    println!("  Best match: {}", classification.best_match);

    // 5. For metagenomic samples, use the Bayesian model for mixtures
    println!("Processing metagenomic sample...");

    // (In a real implementation, we would extract observed k-mer profiles)
    let observed_profile = vec![/* k-mer counts */];

    // Create strain mixture model
    let strain_signatures_matrix = build_signature_matrix(&classifier.references);
    let strain_ids = classifier
        .references
        .iter()
        .map(|r| r.taxon_id.clone())
        .collect();

    let mut mixture_model = StrainMixtureModel::new(
        strain_signatures_matrix,
        strain_ids,
        None, // Use default prior
        None, // Use default iterations
        None, // Use random seed
    )?;

    // Estimate strain abundances
    let result = mixture_model.estimate_abundances(&observed_profile)?;

    println!("Strain abundances:");
    for (strain_id, (abundance, confidence)) in &result.abundances {
        println!(
            "  {}: {:.2}% (±{:.2}%)",
            strain_id,
            abundance * 100.0,
            confidence * 100.0
        );
    }

    Ok(())
}

// Helper function to build signature matrix
fn build_signature_matrix(references: &[MultiResolutionSignature]) -> ndarray::Array2<f64> {
    // This would extract signature features into a matrix
    // Simplified implementation
    ndarray::Array2::zeros((100, references.len()))
}
