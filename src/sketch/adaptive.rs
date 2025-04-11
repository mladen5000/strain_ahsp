//! Adaptive sketching methods.
//!
//! This module could implement techniques like Scaled MinHash or other
//! adaptive methods where the sketch size isn't fixed beforehand but
//! adapts based on the data or a scaling factor. This allows comparing
//! datasets of vastly different sizes more accurately than fixed-size MinHash.

use crate::sketch::signature::Signature;
use crate::sketch::Sketcher; // Implement the common Sketcher trait
use anyhow::Result;
use needletail::parser::FastaReader; // Add this import at the top with other imports
use needletail::parser::SequenceRecord;
use needletail::Sequence;
use std::collections::hash_map::DefaultHasher;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};

/// Structure for creating adaptive sketches (e.g., Scaled MinHash).
#[derive(Debug, Clone)]
pub struct AdaptiveSketcher {
    scaling_factor: u64, // Determines the fraction of hashes to keep (e.g., keep hashes < MAX_HASH / scaling_factor)
    kmer_size: usize,
    // Potentially track max hash value used if needed for specific algorithms
    // max_hash_value: u64,
}

impl AdaptiveSketcher {
    /// Creates a new AdaptiveSketcher.
    ///
    /// # Arguments
    /// * `scaling_factor` - Controls the sketch size. Larger values mean smaller sketches.
    /// * `kmer_size` - The k-mer length.
    pub fn new(scaling_factor: u64, kmer_size: usize) -> Result<Self> {
        if scaling_factor == 0 {
            return Err(anyhow::anyhow!("Scaling factor must be greater than 0."));
        }
        if kmer_size == 0 {
            return Err(anyhow::anyhow!("K-mer size must be greater than 0."));
        }
        Ok(AdaptiveSketcher {
            scaling_factor,
            kmer_size,
            // max_hash_value: u64::MAX / scaling_factor, // Precompute threshold
        })
    }

    /// Calculates a hash value for a k-mer.
    fn hash_kmer(&self, kmer: &[u8]) -> u64 {
        let mut hasher = DefaultHasher::new();
        kmer.hash(&mut hasher);
        hasher.finish()
    }
}

impl Sketcher for AdaptiveSketcher {
    /// Creates an adaptive signature (e.g., Scaled MinHash) for a sequence record.
    /// Keeps all unique k-mer hashes that fall below a threshold determined by the scaling factor.
    fn sketch_sequence(&self, record: &SequenceRecord) -> Result<Signature> {
        // The threshold for keeping hashes
        let threshold = u64::MAX / self.scaling_factor;

        // Note: num_hashes in the signature is less meaningful here, it reflects the *actual*
        // number of hashes kept, which varies. We store 0 or the actual count.
        let mut signature = Signature::new(
            "scaled_minhash".to_string(), // Or other adaptive method name
            self.kmer_size,
            0, // Initial num_hashes is 0, will be updated
        );
        signature.name = Some(String::from_utf8_lossy(record.id()).into_owned());

        let seq = record.sequence();
        let mut kept_hashes = std::collections::HashSet::new(); // Use HashSet to store unique hashes below threshold

        // 1. Generate canonical k-mers and hash them
        // TODO: Use the actual CanonicalKmerIter from bio::kmers when implemented correctly.
        // Manual canonicalization for now:
        for i in 0..=(seq.len().saturating_sub(self.kmer_size)) {
            let kmer = &seq[i..i + self.kmer_size];
            if kmer.iter().any(|&b| !crate::bio::is_valid_base(b)) {
                continue; // Skip k-mers with invalid bases
            }
            let rc = crate::bio::reverse_complement(kmer);
            let canonical_kmer = if kmer < &rc[..] { kmer } else { &rc[..] };

            let hash_value = self.hash_kmer(canonical_kmer);

            // 2. Keep hash if it's below the threshold
            if hash_value < threshold {
                kept_hashes.insert(hash_value);
            }
        }

        // 3. Store the kept hashes in the signature
        signature.hashes = kept_hashes.into_iter().collect();
        signature.hashes.sort_unstable(); // Keep sorted for consistency
        signature.num_hashes = signature.hashes.len(); // Update num_hashes to actual count

        Ok(signature)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use needletail::{parse_fastx_file, parser::SequenceRecord};
    use std::io::Cursor;
    use std::io::{Seek, Write}; // For writing to tempfile in tests

    // Helper to create a SequenceRecord from a string
    fn seq_rec(id: &str, seq: &str) -> SequenceRecord<'static> {
        // Create a temporary file in memory to avoid filesystem dependencies
        let fasta_data = format!(">{}\n{}\n", id, seq);
        let cursor = std::io::Cursor::new(fasta_data.into_bytes());

        // Instead of trying to construct the SequenceRecord manually,
        // we'll use the proper parser but with a memory leak to ensure static lifetime
        let boxed_reader = Box::new(needletail::parse_fastx_reader(cursor).unwrap());

        // Convert the boxed reader to a leaked pointer with 'static lifetime
        // This is unsafe and causes a memory leak, but necessary for the 'static requirement
        let static_reader = Box::leak(boxed_reader);

        // Get the first record from the reader
        match static_reader.next() {
            Some(Ok(record)) => record,
            Some(Err(e)) => panic!("Failed to parse test sequence: {}", e),
            None => panic!("No sequence found in test data"),
        }
    }

    #[test]
    fn test_adaptive_sketcher_new() {
        let sketcher = AdaptiveSketcher::new(1000, 21).unwrap();
        assert_eq!(sketcher.scaling_factor, 1000);
        assert_eq!(sketcher.kmer_size, 21);
    }

    #[test]
    fn test_adaptive_sketcher_new_invalid() {
        assert!(AdaptiveSketcher::new(0, 21).is_err());
        assert!(AdaptiveSketcher::new(1000, 0).is_err());
    }

    #[test]
    fn test_adaptive_sketch_sequence() {
        // This test is difficult without controlling the hash function output.
        // We can check that it runs and produces a signature with some hashes.
        // Use a small scaling factor to increase the chance of keeping hashes.
        let scaling_factor = 10; // Keep hashes roughly < u64::MAX / 10
        let sketcher = AdaptiveSketcher::new(scaling_factor, 3).unwrap();
        let record = seq_rec("test_seq", "ACGTTACGTACGTACTG"); // Longer sequence

        let signature_res = sketcher.sketch_sequence(&record);
        assert!(signature_res.is_ok());
        let signature = signature_res.unwrap();

        assert_eq!(signature.algorithm, "scaled_minhash");
        assert_eq!(signature.kmer_size, 3);
        assert_eq!(signature.name.unwrap(), "test_seq");

        // Check that the number of hashes matches the length of the hash vector
        assert_eq!(signature.num_hashes, signature.hashes.len());
        // Check that *some* hashes were likely kept (probabilistic)
        // This might fail occasionally if all hashes happen to be large
        assert!(
            !signature.hashes.is_empty(),
            "Expected some hashes to be kept with scaling factor {}",
            scaling_factor
        );

        // Verify all kept hashes are below the threshold
        let threshold = u64::MAX / scaling_factor;
        for &hash_val in &signature.hashes {
            assert!(hash_val < threshold);
        }

        // Check that hashes are sorted
        let mut sorted_hashes = signature.hashes.clone();
        sorted_hashes.sort_unstable();
        assert_eq!(signature.hashes, sorted_hashes);
    }

    #[test]
    fn test_adaptive_sketch_empty() {
        let sketcher = AdaptiveSketcher::new(100, 3).unwrap();
        let record = seq_rec("empty", "");
        let signature = sketcher.sketch_sequence(&record).unwrap();
        assert!(signature.hashes.is_empty());
        assert_eq!(signature.num_hashes, 0);
    }

    #[test]
    fn test_adaptive_sketch_k_too_large() {
        let sketcher = AdaptiveSketcher::new(100, 5).unwrap();
        let record = seq_rec("short", "ACGT");
        let signature = sketcher.sketch_sequence(&record).unwrap();
        assert!(signature.hashes.is_empty());
        assert_eq!(signature.num_hashes, 0);
    }
}

/// Classifier using adaptive sketches for taxonomic classification.
/// Implements an algorithm to determine taxonomic identity.
#[derive(Debug)]
pub struct AdaptiveClassifier {
    // Database of reference sketches
    reference_sketches: HashMap<String, Signature>,
    // Scaling factor used for the sketches
    scaling_factor: u64,
    // Minimum similarity threshold to report a match
    min_similarity: f64,
}

impl AdaptiveClassifier {
    /// Creates a new classifier with reference signatures.
    pub fn new(
        reference_sketches: HashMap<String, Signature>,
        scaling_factor: u64,
        min_similarity: f64,
    ) -> Self {
        AdaptiveClassifier {
            reference_sketches,
            scaling_factor,
            min_similarity,
        }
    }

    /// Creates a new classifier with empty references.
    pub fn empty(scaling_factor: u64, min_similarity: f64) -> Self {
        AdaptiveClassifier {
            reference_sketches: HashMap::new(),
            scaling_factor,
            min_similarity,
        }
    }

    /// Adds a reference signature to the classifier.
    pub fn add_reference(&mut self, id: String, signature: Signature) {
        self.reference_sketches.insert(id, signature);
    }

    /// Classifies a query signature against the reference database.
    ///
    /// # Arguments
    ///
    /// * `query_signature` - The signature to classify
    ///
    /// # Returns
    ///
    /// A Vec of (reference ID, similarity score) pairs, sorted by descending similarity
    pub fn classify(&self, query_signature: &Signature) -> Vec<(String, f64)> {
        let mut results = Vec::new();

        for (ref_id, ref_sig) in &self.reference_sketches {
            let similarity = query_signature.jaccard_similarity(ref_sig);

            if similarity >= self.min_similarity {
                results.push((ref_id.clone(), similarity));
            }
        }

        // Sort by similarity (descending)
        results.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

        results
    }

    /// Returns the number of reference sketches in the classifier.
    pub fn reference_count(&self) -> usize {
        self.reference_sketches.len()
    }

    /// Returns the scaling factor used for the sketches.
    pub fn scaling_factor(&self) -> u64 {
        self.scaling_factor
    }

    /// Returns the minimum similarity threshold.
    pub fn min_similarity(&self) -> f64 {
        self.min_similarity
    }
}
