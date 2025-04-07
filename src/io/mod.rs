//! Input/Output operations module.
//!
//! Handles reading data (like FASTQ files, metadata) and writing
//! results (like count tables, analysis outputs).

pub mod fastq; // Sub-module specifically for FASTQ handling

use crate::count_table::CountTable;
use crate::stats::AnalysisResults; // Assuming stats module defines this
use anyhow::Result;
use csv; // Using the csv crate
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Writes analysis results to a CSV file.
///
/// # Arguments
///
/// * `results` - The analysis results structure to write.
/// * `output_path` - The path to the output CSV file.
///
/// # Returns
///
/// * `Result<()>` - Ok(()) if writing was successful, or an error.
pub fn write_results(results: &AnalysisResults, output_path: &str) -> Result<()> {
    let path = Path::new(output_path);
    let file = File::create(path)?;
    let mut writer = csv::Writer::from_writer(BufWriter::new(file));

    // Write header row - Adjust based on AnalysisResults structure
    // Example header:
    writer.write_record(&[
        "feature_id",
        "base_mean", // Average normalized count
        "log2_fold_change",
        "std_error", // Standard error of log2FoldChange
        "stat",      // Wald statistic or similar
        "p_value",
        "p_adjusted", // Adjusted p-value (e.g., Benjamini-Hochberg)
    ])?;

    // Iterate through results and write each row
    // This depends heavily on the structure of AnalysisResults
    for result_item in results.iter() {
        // Assuming results is iterable
        // TODO: Extract data from result_item based on its definition
        let feature_id = &result_item.feature_id;
        let base_mean = result_item.base_mean.to_string();
        let log2fc = result_item
            .log2_fold_change
            .map_or("NA".to_string(), |v| v.to_string());
        let stderr = result_item
            .std_error
            .map_or("NA".to_string(), |v| v.to_string());
        let stat = result_item
            .statistic
            .map_or("NA".to_string(), |v| v.to_string());
        let pval = result_item
            .p_value
            .map_or("NA".to_string(), |v| v.to_string());
        let padj = result_item
            .p_adjusted
            .map_or("NA".to_string(), |v| v.to_string());

        writer.write_record(&[
            feature_id, &base_mean, &log2fc, &stderr, &stat, &pval, &padj,
        ])?;
    }

    writer.flush()?; // Ensure all data is written to the file
    Ok(())
}

/// Writes a CountTable to a CSV file.
///
/// # Arguments
///
/// * `table` - The CountTable to write.
/// * `output_path` - The path to the output CSV file.
///
/// # Returns
///
/// * `Result<()>` - Ok(()) if writing was successful, or an error.
pub fn write_count_table(table: &CountTable, output_path: &str) -> Result<()> {
    let path = Path::new(output_path);
    let file = File::create(path)?;
    let mut writer = csv::Writer::from_writer(BufWriter::new(file));

    // Prepare header: "Feature" followed by sample names
    let mut header = vec!["Feature".to_string()];
    header.extend(table.sample_names().iter().cloned());
    writer.write_record(&header)?;

    // Write rows: feature name followed by counts for each sample
    let counts = table.counts_matrix();
    let (n_features, n_samples) = table.dimensions();
    let feature_names = table.feature_names();

    for r in 0..n_features {
        let mut record = Vec::with_capacity(n_samples + 1);
        record.push(feature_names[r].clone()); // Feature name first
        for c in 0..n_samples {
            record.push(counts[[r, c]].to_string()); // Add count for each sample
        }
        writer.write_record(&record)?;
    }

    writer.flush()?;
    Ok(())
}

/// Reads metadata from a file (typically CSV format).
///
/// # Arguments
///
/// * `metadata_path` - Path to the metadata file.
///
/// # Returns
///
/// * `Result<crate::metadata::Metadata>` - Loaded metadata structure or an error.
pub fn read_metadata(metadata_path: &str) -> Result<crate::metadata::Metadata> {
    crate::metadata::load_metadata(metadata_path)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::count_table::CountTable;
    use crate::stats::DifferentialResult; // Assuming this struct exists
    use ndarray::arr2;
    use std::fs;
    use tempfile::tempdir;

    // Helper to create a simple CountTable for testing
    fn create_test_count_table() -> CountTable {
        let counts = arr2(&[[10.0, 20.0], [5.0, 0.0]]);
        let feature_names = vec!["GeneA", "GeneB"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let sample_names = vec!["Sample1", "Sample2"]
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

    // Helper to create dummy AnalysisResults
    fn create_test_analysis_results() -> AnalysisResults {
        vec![
            DifferentialResult {
                feature_id: "GeneA".to_string(),
                base_mean: 15.0,
                log2_fold_change: Some(1.0),
                std_error: Some(0.5),
                statistic: Some(2.0),
                p_value: Some(0.05),
                p_adjusted: Some(0.1),
            },
            DifferentialResult {
                feature_id: "GeneB".to_string(),
                base_mean: 2.5,
                log2_fold_change: None, // Example with missing values
                std_error: None,
                statistic: None,
                p_value: None,
                p_adjusted: None,
            },
        ]
    }

    #[test]
    fn test_write_count_table_csv() {
        let table = create_test_count_table();
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("counts.csv");
        let output_path_str = file_path.to_str().unwrap();

        write_count_table(&table, output_path_str).unwrap();

        let content = fs::read_to_string(file_path).unwrap();
        let expected_content = "\
Feature,Sample1,Sample2\n\
GeneA,10.0,20.0\n\
GeneB,5.0,0.0\n";
        assert_eq!(content, expected_content);

        dir.close().unwrap();
    }

    #[test]
    fn test_write_results_csv() {
        let results = create_test_analysis_results();
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("results.csv");
        let output_path_str = file_path.to_str().unwrap();

        write_results(&results, output_path_str).unwrap();

        let content = fs::read_to_string(file_path).unwrap();
        let expected_content = "\
feature_id,base_mean,log2_fold_change,std_error,stat,p_value,p_adjusted\n\
GeneA,15.0,1.0,0.5,2.0,0.05,0.1\n\
GeneB,2.5,NA,NA,NA,NA,NA\n"; // Note NA for None values
        assert_eq!(content, expected_content);

        dir.close().unwrap();
    }
}
