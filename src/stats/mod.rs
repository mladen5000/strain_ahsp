//! Statistical analysis module.
//!
//! This module contains functions for performing statistical tests,
//! likely focusing on differential abundance analysis similar to DESeq2.
//! It might also include general statistical utilities or other analysis types.

pub mod bayesian; // Sub-module for Bayesian statistical methods
pub mod deconvolution;

pub use bayesian::StrainMixtureModel;
pub use deconvolution::StrainDeconvolution;

use crate::count_table::CountTable;
use crate::metadata::load_metadata;
use anyhow::Result;
use serde::{Deserialize, Serialize};

/// Represents the results of a differential abundance analysis for a single feature.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DifferentialResult {
    pub feature_id: String,
    pub base_mean: f64,                // Mean normalized count across all samples
    pub log2_fold_change: Option<f64>, // Estimated effect size
    pub std_error: Option<f64>,        // Standard error of log2FoldChange estimate
    pub statistic: Option<f64>,        // Wald statistic or similar test statistic
    pub p_value: Option<f64>,          // Raw p-value from the test
    pub p_adjusted: Option<f64>,       // Adjusted p-value (e.g., Benjamini-Hochberg)
}

/// Type alias for the collection of results from an analysis.
pub type AnalysisResults = Vec<DifferentialResult>;

/// Re-export Metadata from metadata module for backward compatibility
pub use crate::metadata::Metadata; // metadata::Metadata as SampleMetadata

/// Runs the main differential abundance analysis (e.g., DESeq2-like).
///
/// # Arguments
///
/// * `normalized_table` - The CountTable with normalized counts.
/// * `metadata_path` - Path to the metadata file describing samples and conditions.
///
/// # Returns
///
/// * `Result<AnalysisResults>` - A vector of results for each feature, or an error.
pub fn run_deseq2_like_analysis(
    normalized_table: &CountTable,
    metadata_path: &Option<String>,
) -> Result<AnalysisResults> {
    // 1. Load and validate metadata
    let metadata = match metadata_path {
        Some(path) => load_metadata(path)?,
        None => {
            return Err(anyhow::anyhow!(
                "Metadata file is required for differential analysis."
            ))
        }
    };
    validate_metadata(normalized_table, &metadata)?;

    // 2. TODO: Implement the core DESeq2-like algorithm. This is complex and involves:
    //    a. Estimating size factors (often done during normalization, but might be re-checked).
    //    b. Estimating dispersion for each feature (variance estimation, crucial step).
    //       - Fit dispersion trend (mean-dispersion relationship).
    //       - Shrink feature-wise estimates towards the trend (empirical Bayes shrinkage).
    //    c. Fitting a Generalized Linear Model (GLM) for each feature.
    //       - Typically Negative Binomial GLM: count ~ condition + other_covariates.
    //    d. Performing hypothesis testing on model coefficients (e.g., Wald test for condition effect).
    //    e. Adjusting p-values for multiple testing (e.g., Benjamini-Hochberg).

    unimplemented!("Core DESeq2-like analysis (dispersion estimation, GLM fitting, testing) needs implementation.");

    // Placeholder return
    // Ok(Vec::new())
}

/// Loads metadata from a file (e.g., CSV).
///
/// # Arguments
/// * `path` - Path to the metadata file.
///
/// # Returns
/// * `Result<Metadata>` - Loaded metadata or an error.
/// IDK
fn load_some_metadata(path: &str) -> Result<Metadata> {
    load_metadata(path)
}

/// Validates that the metadata matches the samples in the CountTable.
///
/// # Arguments
/// * `table` - The CountTable.
/// * `metadata` - The loaded Metadata.
///
/// # Returns
/// * `Result<()>` - Ok(()) if valid, or an error describing mismatches.
fn validate_metadata(table: &CountTable, metadata: &Metadata) -> Result<()> {
    let table_samples: std::collections::HashSet<_> =
        table.sample_names().iter().cloned().collect();
    let metadata_samples: std::collections::HashSet<_> =
        metadata.sample_info.keys().cloned().collect();

    let table_only: Vec<_> = table_samples.difference(&metadata_samples).collect();
    let metadata_only: Vec<_> = metadata_samples.difference(&table_samples).collect();

    let mut errors = Vec::new();
    if !table_only.is_empty() {
        errors.push(format!(
            "Samples found in count table but not in metadata: {:?}",
            table_only
        ));
    }
    if !metadata_only.is_empty() {
        errors.push(format!(
            "Samples found in metadata but not in count table: {:?}",
            metadata_only
        ));
    }

    // Check for samples with empty conditions in metadata that are present in the table
    for sample_name in table.sample_names() {
        if let Some(condition) = metadata.condition_map.get(sample_name) {
            if condition.is_empty() {
                errors.push(format!(
                    "Sample '{}' has an empty condition in the metadata.",
                    sample_name
                ));
            }
        }
        // If sample not in metadata, it's already covered by metadata_only check
    }

    if !errors.is_empty() {
        Err(anyhow::anyhow!(
            "Metadata validation failed:\n- {}",
            errors.join("\n- ")
        ))
    } else if table_samples.is_empty() {
        Err(anyhow::anyhow!("Count table contains no samples."))
    } else {
        log::info!("Metadata validated successfully against count table samples.");
        Ok(())
    }
}

/// Adjusts p-values for multiple testing using Benjamini-Hochberg method.
///
/// # Arguments
/// * `results` - A mutable slice of DifferentialResult structs containing raw p-values.
pub fn adjust_pvalues_bh(results: &mut [DifferentialResult]) {
    // Sort results by p-value, keeping track of original indices
    let mut indexed_results: Vec<(usize, Option<f64>)> = results
        .iter()
        .enumerate()
        .map(|(i, r)| (i, r.p_value))
        .collect();

    // Sort by p-value, putting None (NA) values last
    indexed_results.sort_unstable_by(|a, b| match (a.1, b.1) {
        (Some(pa), Some(pb)) => pa.partial_cmp(&pb).unwrap_or(std::cmp::Ordering::Equal),
        (Some(_), None) => std::cmp::Ordering::Less,
        (None, Some(_)) => std::cmp::Ordering::Greater,
        (None, None) => std::cmp::Ordering::Equal,
    });

    let m = indexed_results.iter().filter(|(_, p)| p.is_some()).count(); // Number of tests performed (non-NA p-values)
    let mut last_padj = 1.0; // Start from the highest possible adjusted p-value

    // Iterate downwards through sorted p-values
    for (rank, (original_index, p_value_opt)) in indexed_results.iter().enumerate().rev() {
        if let Some(p_value) = p_value_opt {
            if m == 0 {
                // Avoid division by zero if no tests were done
                results[*original_index].p_adjusted = Some(1.0);
                continue;
            }
            // Calculate BH adjusted p-value: p * m / rank
            // Rank here is 1-based index of the sorted p-value
            let rank_1_based = rank + 1;
            let padj = (p_value * m as f64) / rank_1_based as f64;

            // Enforce monotonicity: adjusted p-value cannot be greater than the next highest
            let current_padj = padj.min(last_padj).min(1.0); // Ensure it doesn't exceed 1.0

            results[*original_index].p_adjusted = Some(current_padj);
            last_padj = current_padj; // Update the last adjusted p-value seen
        } else {
            // If p-value was None, adjusted p-value is also None
            results[*original_index].p_adjusted = None;
            // Keep last_padj as is, effectively carrying forward the previous value for monotonicity check
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::count_table::CountTable;
    use ndarray::arr2;
    use std::fs::File;
    use std::io::Write;
    use tempfile::tempdir;

    // Helper to create a simple CountTable for testing metadata validation
    fn create_meta_test_table() -> CountTable {
        let counts = arr2(&[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
        let feature_names: Vec<String> = vec!["F1", "F2"].iter().map(|s| s.to_string()).collect();
        let sample_names: Vec<String> = vec!["S1", "S2", "S3"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let feature_map = feature_names
            .iter()
            .enumerate()
            .map(|(i, n)| (n.clone(), i))
            .collect();
        let sample_map = sample_names
            .iter()
            .enumerate()
            .map(|(i, n)| (n.clone(), i))
            .collect();
        CountTable {
            counts,
            feature_names,
            feature_map,
            sample_names,
            sample_map,
        }
    }

    // Helper to create a dummy metadata file
    fn create_dummy_metadata(path: &std::path::Path, content: &str) {
        let mut file = File::create(path).unwrap();
        writeln!(file, "{}", content).unwrap();
    }

    #[test]
    fn test_load_metadata_ok() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("meta.csv");
        create_dummy_metadata(
            &file_path,
            "SampleID,Condition\nS1,Control\nS2,Treatment\nS3,Control",
        );

        let metadata_res = load_metadata(file_path.to_str().unwrap());
        assert!(metadata_res.is_ok());

        let metadata = metadata_res.unwrap();
        assert_eq!(metadata.condition_map.len(), 3);
        assert_eq!(
            metadata.condition_map.get("S1"),
            Some(&"Control".to_string())
        );
        assert_eq!(
            metadata.condition_map.get("S2"),
            Some(&"Treatment".to_string())
        );
        assert_eq!(
            metadata.condition_map.get("S3"),
            Some(&"Control".to_string())
        );
    }
}
