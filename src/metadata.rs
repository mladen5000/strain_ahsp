//! Metadata handling module.
//!
//! This module provides structures and functions for working with sample metadata,
//! including experimental design and sample information.

use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;

/// Represents metadata for a collection of samples.
/// Contains information about sample conditions, groups, and other relevant metadata.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Metadata {
    /// Maps sample IDs to their conditions/groups
    pub condition_map: HashMap<String, String>,
    /// Additional sample metadata, extensible for different experiment designs
    pub sample_attributes: HashMap<String, HashMap<String, String>>,
}

impl Metadata {
    /// Creates a new, empty Metadata structure
    pub fn new() -> Self {
        Metadata {
            condition_map: HashMap::new(),
            sample_attributes: HashMap::new(),
        }
    }

    /// Adds a sample with its condition to the metadata
    pub fn add_sample(&mut self, sample_id: &str, condition: &str) {
        self.condition_map
            .insert(sample_id.to_string(), condition.to_string());
    }

    /// Adds an attribute for a specific sample
    pub fn add_sample_attribute(&mut self, sample_id: &str, attribute: &str, value: &str) {
        let sample_map = self
            .sample_attributes
            .entry(sample_id.to_string())
            .or_insert_with(HashMap::new);
        sample_map.insert(attribute.to_string(), value.to_string());
    }

    /// Returns all unique conditions present in the metadata
    pub fn get_conditions(&self) -> Vec<String> {
        let mut conditions: Vec<String> = self.condition_map.values().cloned().collect();

        // Remove duplicates
        conditions.sort();
        conditions.dedup();

        conditions
    }

    /// Returns the number of samples in the metadata
    pub fn sample_count(&self) -> usize {
        self.condition_map.len()
    }
}

/// Loads metadata from a CSV file.
///
/// # Arguments
///
/// * `path` - Path to the metadata CSV file
///
/// # Returns
///
/// * `Result<Metadata>` - Metadata structure or error
pub fn load_metadata(path: &str) -> Result<Metadata> {
    let path = Path::new(path);
    let mut rdr = csv::Reader::from_path(path)?;
    let mut metadata = Metadata::new();

    // Find column indices for required fields
    let headers = rdr.headers()?.clone();
    let sample_col = headers.iter().position(|h| {
        h.trim().eq_ignore_ascii_case("sampleid") || h.trim().eq_ignore_ascii_case("sample")
    });
    let condition_col = headers.iter().position(|h| {
        h.trim().eq_ignore_ascii_case("condition") || h.trim().eq_ignore_ascii_case("group")
    });

    let sample_col = sample_col
        .ok_or_else(|| anyhow::anyhow!("Metadata CSV missing 'SampleID'/'Sample' column"))?;
    let condition_col = condition_col
        .ok_or_else(|| anyhow::anyhow!("Metadata CSV missing 'Condition'/'Group' column"))?;

    // Process each row
    for result in rdr.records() {
        let record = result?;
        let sample_id = record
            .get(sample_col)
            .ok_or_else(|| anyhow::anyhow!("Missing sample ID in metadata row"))?
            .trim()
            .to_string();
        let condition = record
            .get(condition_col)
            .ok_or_else(|| anyhow::anyhow!("Missing condition in metadata row"))?
            .trim()
            .to_string();

        if sample_id.is_empty() {
            log::warn!("Skipping metadata row with empty sample ID.");
            continue;
        }
        if condition.is_empty() {
            log::warn!("Sample '{}' has an empty condition in metadata.", sample_id);
            // Store it anyway, validation will catch it later
        }

        metadata.add_sample(&sample_id, &condition);

        // Process remaining columns as additional attributes
        for (i, field) in record.iter().enumerate() {
            // Skip sample ID and condition columns which we've already processed
            if i != sample_col && i != condition_col {
                let header = headers.get(i).unwrap_or(&"unknown");
                metadata.add_sample_attribute(&sample_id, header, field);
            }
        }
    }

    if metadata.sample_count() == 0 {
        return Err(anyhow::anyhow!(
            "No valid sample entries found in metadata file '{}'",
            path.display()
        ));
    }

    Ok(metadata)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Write;
    use tempfile::tempdir;

    fn create_test_metadata_file(path: &std::path::Path, content: &str) {
        let mut file = File::create(path).unwrap();
        writeln!(file, "{}", content).unwrap();
    }

    #[test]
    fn test_load_metadata_basic() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("metadata.csv");
        create_test_metadata_file(
            &file_path,
            "SampleID,Condition,Batch\nS1,Control,B1\nS2,Treatment,B1\nS3,Control,B2",
        );

        let metadata = load_metadata(file_path.to_str().unwrap()).unwrap();

        assert_eq!(metadata.sample_count(), 3);
        assert_eq!(
            metadata.condition_map.get("S1"),
            Some(&"Control".to_string())
        );
        assert_eq!(
            metadata.condition_map.get("S2"),
            Some(&"Treatment".to_string())
        );

        // Check additional attributes
        assert_eq!(
            metadata.sample_attributes.get("S1").unwrap().get("Batch"),
            Some(&"B1".to_string())
        );

        // Check conditions list
        let conditions = metadata.get_conditions();
        assert_eq!(conditions.len(), 2);
        assert!(conditions.contains(&"Control".to_string()));
        assert!(conditions.contains(&"Treatment".to_string()));
    }

    #[test]
    fn test_load_metadata_missing_columns() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("invalid.csv");
        create_test_metadata_file(&file_path, "Sample,OtherField\nS1,Value1\n");

        let result = load_metadata(file_path.to_str().unwrap());
        assert!(result.is_err());

        // Different column names
        let file_path2 = dir.path().join("valid_alt_names.csv");
        create_test_metadata_file(&file_path2, "sample,group\nS1,Control\n");

        let result2 = load_metadata(file_path2.to_str().unwrap());
        assert!(result2.is_ok());
    }
}
