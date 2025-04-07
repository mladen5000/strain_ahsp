//! Defines structures and functions for handling count data.
//!
//! This module likely deals with representing and manipulating tables
//! where rows might be features (genes, k-mers, taxa) and columns
//! are samples.

use anyhow::Result;
use ndarray::{Array, Array2, Axis}; // Using ndarray for matrix operations
use serde::{Deserialize, Serialize};
use std::collections::HashMap; // Or indexmap::IndexMap for ordered keys // For potential serialization

/// Represents a count table.
///
/// Stores counts (e.g., u32 or f64 after normalization) along with
/// mappings for feature names and sample names.
#[derive(Debug, Serialize, Deserialize)]
pub struct CountTable {
    /// The core count data matrix (features x samples).
    pub counts: Array2<f64>, // Use f64 to accommodate normalized values

    /// Mapping from feature index (row) to feature name (e.g., k-mer string, gene ID).
    pub feature_names: Vec<String>,
    pub feature_map: HashMap<String, usize>, // For quick lookup

    /// Mapping from sample index (column) to sample name.
    pub sample_names: Vec<String>,
    pub sample_map: HashMap<String, usize>, // For quick lookup
}

impl CountTable {
    /// Creates a new, empty CountTable.
    pub fn new() -> Self {
        // TODO: Consider initializing with expected dimensions if known.
        CountTable {
            counts: Array2::zeros((0, 0)),
            feature_names: Vec::new(),
            feature_map: HashMap::new(),
            sample_names: Vec::new(),
            sample_map: HashMap::new(),
        }
    }

    /// Builds a CountTable from processed data (e.g., k-mer counts per sample).
    ///
    /// # Arguments
    ///
    /// * `data` - A structure representing the counts per sample. This needs definition.
    ///   For example, it could be `HashMap<String, HashMap<String, u32>>` where
    ///   outer key is sample name, inner key is feature name, value is count.
    ///
    /// # Returns
    ///
    /// * `Result<Self>` - The constructed CountTable or an error.
    pub fn build_from_data(/* data: AppropriateDataStructure */) -> Result<Self> {
        // TODO: Implement the logic to:
        // 1. Collect all unique feature names and sample names.
        // 2. Create the feature_names, feature_map, sample_names, sample_map.
        // 3. Initialize the ndarray::Array2 with zeros based on dimensions.
        // 4. Populate the Array2 with counts from the input data structure.
        // 5. Handle potential errors (e.g., inconsistent data).
        unimplemented!("CountTable::build_from_data needs implementation");
    }

    /// Adds a sample column to the table.
    ///
    /// # Arguments
    ///
    /// * `sample_name` - The name of the sample to add.
    /// * `sample_counts` - A map or vector of counts for this sample (feature -> count).
    pub fn add_sample(
        &mut self,
        sample_name: &str, /* sample_counts: AppropriateDataStructure */
    ) -> Result<()> {
        // TODO: Implement logic to add a new column.
        // - Check if sample already exists.
        // - Add new features if encountered.
        // - Resize the array if necessary.
        // - Fill in the counts for the new sample.
        unimplemented!("CountTable::add_sample needs implementation");
    }

    /// Retrieves the counts for a specific feature.
    pub fn get_feature_counts(&self, feature_name: &str) -> Option<ndarray::ArrayView1<f64>> {
        // TODO: Implement lookup using feature_map and return a view of the row.
        self.feature_map
            .get(feature_name)
            .map(|&idx| self.counts.row(idx))
    }

    /// Retrieves the counts for a specific sample.
    pub fn get_sample_counts(&self, sample_name: &str) -> Option<ndarray::ArrayView1<f64>> {
        // TODO: Implement lookup using sample_map and return a view of the column.
        self.sample_map
            .get(sample_name)
            .map(|&idx| self.counts.column(idx))
    }

    /// Returns the dimensions of the count table (features, samples).
    pub fn dimensions(&self) -> (usize, usize) {
        self.counts.dim()
    }

    /// Returns a reference to the underlying count matrix.
    pub fn counts_matrix(&self) -> &Array2<f64> {
        &self.counts
    }

    /// Returns a mutable reference to the underlying count matrix.
    pub fn counts_matrix_mut(&mut self) -> &mut Array2<f64> {
        &mut self.counts
    }

    /// Returns the list of feature names.
    pub fn feature_names(&self) -> &Vec<String> {
        &self.feature_names
    }

    /// Returns the list of sample names.
    pub fn sample_names(&self) -> &Vec<String> {
        &self.sample_names
    }

    // TODO: Add methods for filtering (by count, by feature presence), merging tables, etc.
    // TODO: Add methods for writing the table to a file (e.g., CSV).
}

// Default implementation for creating an empty table.
impl Default for CountTable {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::arr2;

    #[test]
    fn test_new_count_table() {
        let table = CountTable::new();
        assert_eq!(table.dimensions(), (0, 0));
        assert!(table.feature_names.is_empty());
        assert!(table.sample_names.is_empty());
    }

    // TODO: Add more comprehensive tests once build_from_data and other methods are implemented.
    // Example structure for a test:
    // #[test]
    // fn test_build_simple_table() {
    //     // 1. Create mock input data (e.g., HashMap)
    //     // 2. Call CountTable::build_from_data
    //     // 3. Assert dimensions are correct
    //     // 4. Assert feature/sample names are correct
    //     // 5. Assert specific count values are correct
    //     // assert!(false, "Test not implemented");
    // }
}
