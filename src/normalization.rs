//! Provides functions for normalizing count data.
//!
//! Normalization aims to adjust raw counts to account for differences
//! in sequencing depth, library size, or other factors, making counts
//! comparable across samples.

use crate::count_table::CountTable;
use anyhow::{anyhow, Result};
use log::warn;
use ndarray::{s, Array1, ArrayView1, Axis};
use statrs::statistics::Data; // For median calculation
use statrs::statistics::OrderStatistics; // For median calculation

/// Normalizes the counts in a CountTable using a specified method.
///
/// # Arguments
///
/// * `table` - A mutable reference to the CountTable to normalize.
/// * `method` - A string slice specifying the normalization method
///              (e.g., "median-of-ratios", "tpm", "cpm", "none").
///
/// # Returns
///
/// * `Result<()>` - Ok if normalization was successful, otherwise an error.
pub fn normalize(table: &mut CountTable, method: &str) -> Result<()> {
    match method.to_lowercase().as_str() {
        "median-of-ratios" | "deseq2" => normalize_median_of_ratios(table),
        "tpm" => normalize_tpm(table), // Requires gene lengths - needs modification
        "cpm" => normalize_cpm(table),
        "none" => {
            warn!("No normalization applied.");
            Ok(())
        }
        _ => Err(anyhow!("Unsupported normalization method: {}", method)),
    }
}

/// Normalizes counts using the Median-of-Ratios method (similar to DESeq2).
///
/// 1. Calculate a pseudo-reference sample (geometric mean of counts across samples for each feature).
/// 2. For each sample, calculate the ratio of its counts to the pseudo-reference for each feature.
/// 3. Calculate the median of these ratios for each sample (this is the size factor).
/// 4. Divide the counts in each sample by its size factor.
///
/// # Arguments
///
/// * `table` - A mutable reference to the CountTable.
fn normalize_median_of_ratios(table: &mut CountTable) -> Result<()> {
    let counts = table.counts_matrix();
    let (n_features, n_samples) = counts.dim();
    let sample_names = table.sample_names().to_vec(); // Store sample names upfront

    if n_features == 0 || n_samples == 0 {
        warn!("Count table is empty, skipping median-of-ratios normalization.");
        return Ok(());
    }

    // Calculate geometric mean for each feature across samples, ignoring zeros
    let mut log_counts_sum = Array1::<f64>::zeros(n_features);
    let mut non_zero_counts = Array1::<f64>::zeros(n_features);

    for r in 0..n_features {
        for c in 0..n_samples {
            let count = counts[[r, c]];
            if count > 0.0 {
                log_counts_sum[r] += count.ln();
                non_zero_counts[r] += 1.0;
            }
        }
    }

    let pseudo_reference = log_counts_sum
        .iter()
        .zip(non_zero_counts.iter())
        .map(|(&sum, &count)| {
            if count > 0.0 {
                (sum / count).exp()
            } else {
                0.0
            }
        })
        .collect::<Array1<f64>>();

    // Calculate size factors for each sample
    let mut size_factors = Array1::<f64>::zeros(n_samples);
    for c in 0..n_samples {
        let mut ratios = Vec::new();
        for r in 0..n_features {
            let count = counts[[r, c]];
            let ref_val = pseudo_reference[r];
            // Only consider features with non-zero counts in both sample and reference
            if count > 0.0 && ref_val > 0.0 {
                ratios.push(count / ref_val);
            }
        }

        if ratios.is_empty() {
            warn!("Sample {} has no features with positive counts common with the pseudo-reference. Setting size factor to 1.0.", sample_names[c]);
            size_factors[c] = 1.0; // Or handle as error? Or use total count?
        } else {
            // Calculate median of ratios
            let mut data = Data::new(ratios);
            size_factors[c] = data.median();
            if size_factors[c] <= 0.0 || !size_factors[c].is_finite() {
                warn!("Calculated non-positive or non-finite size factor ({}) for sample {}. Setting to 1.0.", size_factors[c], sample_names[c]);
                size_factors[c] = 1.0; // Fallback if median is zero or invalid
            }
        }
    }

    // Normalize counts by dividing each sample's counts by its size factor
    let mut normalized_counts = table.counts_matrix_mut();
    for c in 0..n_samples {
        let sf = size_factors[c];
        if sf > 0.0 && sf.is_finite() {
            // Ensure size factor is valid
            let mut sample_col = normalized_counts.column_mut(c);
            sample_col /= sf;
        } else {
            warn!(
                "Skipping normalization for sample {} due to invalid size factor {}.",
                sample_names[c], sf
            );
        }
    }

    Ok(())
}

/// Normalizes counts to Counts Per Million (CPM).
/// CPM = (count / total_counts_in_sample) * 1,000,000
///
/// # Arguments
///
/// * `table` - A mutable reference to the CountTable.
fn normalize_cpm(table: &mut CountTable) -> Result<()> {
    let mut counts = table.counts_matrix_mut();
    let library_sizes = counts.sum_axis(Axis(0)); // Sum counts per sample (column)

    if library_sizes.iter().any(|&sum| sum <= 0.0) {
        warn!("Some samples have zero or negative total counts. CPM normalization might produce NaNs or Infs.");
    }

    counts
        .axis_iter_mut(Axis(1)) // Iterate over columns (samples)
        .zip(library_sizes.iter())
        .for_each(|(mut col, &total_counts)| {
            if total_counts > 0.0 {
                col *= 1_000_000.0 / total_counts;
            } else {
                // Handle samples with zero total counts - set to zero? Or leave as is?
                col.fill(0.0);
            }
        });

    Ok(())
}

/// Normalizes counts to Transcripts Per Million (TPM).
/// TPM requires feature lengths (e.g., gene lengths).
/// 1. Divide each count by the length of its feature (in kilobases).
/// 2. Calculate the sum of these length-normalized counts for each sample ("per million" scaling factor).
/// 3. Divide the length-normalized counts by the scaling factor and multiply by 1,000,000.
///
/// # Arguments
///
/// * `table` - A mutable reference to the CountTable.
/// * `feature_lengths` - A slice or map providing the length for each feature. **This needs to be passed in.**
fn normalize_tpm(
    table: &mut CountTable, /*, feature_lengths: &[f64] or HashMap<String, f64> */
) -> Result<()> {
    // TODO: Implement TPM normalization.
    // This requires feature lengths, which are not currently part of CountTable.
    // The function signature needs to be updated to accept lengths.
    // 1. Get feature lengths matching the order in table.feature_names.
    // 2. Divide counts by length in kb (count / (length / 1000.0)).
    // 3. Sum these values per sample.
    // 4. Divide length-normalized counts by the sum and multiply by 1e6.
    Err(anyhow!(
        "TPM normalization requires feature lengths and is not yet implemented."
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::count_table::CountTable;
    use approx::assert_relative_eq;
    use ndarray::{arr2, Array2, Axis}; // For float comparisons

    // Helper to create a simple CountTable for testing
    fn create_test_table() -> CountTable {
        let counts = arr2(&[
            [10.0, 20.0, 30.0], // Feature 1
            [5.0, 0.0, 15.0],   // Feature 2 (zero in sample 2)
            [0.0, 40.0, 60.0],  // Feature 3 (zero in sample 1)
            [2.0, 4.0, 6.0],    // Feature 4
        ]);
        let feature_names: Vec<String> = vec!["F1", "F2", "F3", "F4"]
            .iter()
            .map(|s| s.to_string())
            .collect();
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

    #[test]
    fn test_normalize_cpm() {
        let mut table = create_test_table();
        normalize_cpm(&mut table).unwrap();

        let expected_totals = [17.0, 64.0, 111.0]; // Original totals
        let mut expected_cpm = table.counts_matrix().to_owned();
        for c in 0..3 {
            let total = expected_totals[c];
            if total > 0.0 {
                let mut col = expected_cpm.column_mut(c);
                col *= 1_000_000.0 / total;
            } else {
                expected_cpm.column_mut(c).fill(0.0);
            }
        }

        // Compare arrays element-wise
        let actual = table.counts_matrix();
        for i in 0..actual.nrows() {
            for j in 0..actual.ncols() {
                assert_relative_eq!(actual[[i, j]], expected_cpm[[i, j]], epsilon = 1e-6);
            }
        }
    }

    #[test]
    fn test_normalize_median_of_ratios() {
        // Note: This is a simplified test. Real DESeq2 involves more nuances.
        let mut table = create_test_table();
        // Original:
        // [[10.0, 20.0, 30.0],
        //  [ 5.0,  0.0, 15.0],
        //  [ 0.0, 40.0, 60.0],
        //  [ 2.0,  4.0,  6.0]]

        normalize_median_of_ratios(&mut table).unwrap();

        // Expected pseudo-reference (geometric mean, approx):
        // F1: (10*20*30)^(1/3) ~= 18.17
        // F2: (5*15)^(1/2) ~= 8.66  (ignore zero)
        // F3: (40*60)^(1/2) ~= 48.99 (ignore zero)
        // F4: (2*4*6)^(1/3) ~= 3.63

        // Expected ratios (approx):
        // S1: [10/18.17, 5/8.66, NA, 2/3.63] = [0.55, 0.58, NA, 0.55] -> median ~0.55
        // S2: [20/18.17, NA, 40/48.99, 4/3.63] = [1.10, NA, 0.82, 1.10] -> median ~1.10
        // S3: [30/18.17, 15/8.66, 60/48.99, 6/3.63] = [1.65, 1.73, 1.22, 1.65] -> median ~1.65

        // Expected size factors (approx): [0.55, 1.10, 1.65]

        // Expected normalized counts (original / size_factor, approx):
        let expected = arr2(&[
            [10.0 / 0.55, 20.0 / 1.10, 30.0 / 1.65], // ~ [18.18, 18.18, 18.18]
            [5.0 / 0.55, 0.0 / 1.10, 15.0 / 1.65],   // ~ [ 9.09,  0.00,  9.09]
            [0.0 / 0.55, 40.0 / 1.10, 60.0 / 1.65],  // ~ [ 0.00, 36.36, 36.36]
            [2.0 / 0.55, 4.0 / 1.10, 6.0 / 1.65],    // ~ [ 3.64,  3.64,  3.64]
        ]);

        println!("Normalized Table:\n{:?}", table.counts_matrix());
        println!("Expected Table:\n{:?}", expected);

        // Compare arrays element-wise
        let actual = table.counts_matrix();
        for i in 0..actual.nrows() {
            for j in 0..actual.ncols() {
                assert_relative_eq!(actual[[i, j]], expected[[i, j]], epsilon = 1e-2);
            }
        }
    }

    #[test]
    fn test_normalize_none() {
        let mut table = create_test_table();
        let original_counts = table.counts_matrix().to_owned();
        normalize(&mut table, "none").unwrap();
        assert_eq!(table.counts_matrix(), &original_counts);
    }

    #[test]
    fn test_unsupported_method() {
        let mut table = create_test_table();
        let result = normalize(&mut table, "unknown_method");
        assert!(result.is_err());
    }

    #[test]
    fn test_tpm_unimplemented() {
        let mut table = create_test_table();
        let result = normalize_tpm(&mut table); // No lengths provided
        assert!(result.is_err());
    }
}
