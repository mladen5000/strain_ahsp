//! Defines structures and functions related to sequence signatures.
//!
//! This module is used in conjunction with sketching techniques
//! (like MinHash) to represent sequences or datasets compactly.
//! 
//! Note: Some functionality is now re-exported from sketch/signature module
//! to maintain API compatibility.

use anyhow::Result;
use serde::{Deserialize, Serialize};

// Assuming MinHash is used (based on sketch/minhash.rs)
// A signature often consists of a collection of hash values.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Signature {
    /// The algorithm used (e.g., "minhash").
    pub algorithm: String,
    /// The k-mer size used for generating the signature.
    pub kmer_size: usize,
    /// The number of hash values (or maximum size for adaptive methods).
    pub num_hashes: usize, // Or max_size for adaptive
    /// The actual hash values comprising the signature.
    pub hashes: Vec<u64>, // Assuming u64 hashes, adjust if different
    /// Optional: Name or identifier for the sequence/dataset.
    pub name: Option<String>,
    /// Optional: Filename from which the signature was derived.
    pub filename: Option<String>,
    // Add other metadata as needed (e.g., molecule type DNA/protein, alphabet)
}

impl Signature {
    /// Creates a new, empty signature placeholder.
    pub fn new(algorithm: String, kmer_size: usize, num_hashes: usize) -> Self {
        Signature {
            algorithm,
            kmer_size,
            num_hashes,
            hashes: Vec::with_capacity(num_hashes), // Pre-allocate if fixed size
            name: None,
            filename: None,
        }
    }

    /// Adds a hash value to the signature.
    /// For fixed-size MinHash, this might involve checks or specific insertion logic.
    /// For adaptive methods, it might just append.
    pub fn add_hash(&mut self, hash: u64) {
        // TODO: Implement logic based on the specific sketching algorithm.
        // For simple MinHash, ensure we don't exceed num_hashes and potentially keep sorted.
        if self.hashes.len() < self.num_hashes {
            self.hashes.push(hash);
            // Optionally keep sorted for faster comparisons later
            // self.hashes.sort_unstable();
        } else {
            // If using standard MinHash (fixed size), we might need to replace
            // the largest hash if the new one is smaller.
            // This requires keeping the hashes sorted or finding the max.
            // For now, just ignoring overflow for simplicity.
            log::warn!(
                "Signature hash capacity ({}) reached. Ignoring new hash {}.",
                self.num_hashes,
                hash
            );
        }
    }

    /// Merges another signature into this one.
    /// Typically used in MinHash to combine sketches.
    pub fn merge(&mut self, other: &Signature) -> Result<()> {
        // TODO: Implement signature merging logic.
        // - Check compatibility (algorithm, k-mer size).
        // - Combine hash sets (e.g., take the union and keep the smallest `num_hashes`).
        if self.algorithm != other.algorithm || self.kmer_size != other.kmer_size {
            return Err(anyhow::anyhow!(
                "Cannot merge signatures with different algorithms or k-mer sizes."
            ));
        }
        unimplemented!("Signature::merge needs implementation.");
    }

    /// Calculates the Jaccard similarity between this signature and another.
    /// Jaccard = |Intersection| / |Union|
    pub fn jaccard(&self, other: &Signature) -> Result<f64> {
        if self.algorithm != other.algorithm || self.kmer_size != other.kmer_size {
            return Err(anyhow::anyhow!(
                "Cannot compare signatures with different algorithms or k-mer sizes."
            ));
        }
        if self.num_hashes != other.num_hashes {
            // Jaccard estimation with standard MinHash assumes equal sketch sizes.
            // Comparison might still be possible but interpretation differs.
            log::warn!("Comparing signatures with different numbers of hashes ({} vs {}). Jaccard estimate might be less accurate.", self.num_hashes, other.num_hashes);
        }

        if self.hashes.is_empty() && other.hashes.is_empty() {
            return Ok(1.0); // Jaccard of two empty sets is often defined as 1
        }
        if self.hashes.is_empty() || other.hashes.is_empty() {
            return Ok(0.0); // Jaccard with an empty set is 0
        }

        // TODO: Implement efficient intersection calculation.
        // Assumes hashes are unique within each signature (typical for MinHash).
        // Using a HashSet for intersection is simple but might not be the most performant.
        use std::collections::HashSet;
        let set1: HashSet<_> = self.hashes.iter().collect();
        let set2: HashSet<_> = other.hashes.iter().collect();

        let intersection_size = set1.intersection(&set2).count();
        let union_size = set1.union(&set2).count(); // Or set1.len() + set2.len() - intersection_size

        if union_size == 0 {
            Ok(1.0) // Avoid division by zero if somehow union is empty (shouldn't happen if intersection works)
        } else {
            Ok(intersection_size as f64 / union_size as f64)
        }
    }

    // TODO: Add methods for saving/loading signatures (e.g., to/from JSON or binary format).
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;
    use approx::assert_relative_eq;

    fn create_test_sig(name: &str, hashes: Vec<u64>) -> Signature {
        Signature {
            algorithm: "minhash".to_string(),
            kmer_size: 21,
            num_hashes: hashes.len(), // Use actual length for testing simplicity
            hashes,
            name: Some(name.to_string()),
            filename: None,
        }
    }

    #[test]
    fn test_signature_new() {
        let sig = Signature::new("minhash".to_string(), 21, 100);
        assert_eq!(sig.algorithm, "minhash");
        assert_eq!(sig.kmer_size, 21);
        assert_eq!(sig.num_hashes, 100);
        assert!(sig.hashes.is_empty());
    }

    #[test]
    fn test_signature_add_hash() {
        let mut sig = Signature::new("minhash".to_string(), 21, 2);
        sig.add_hash(100);
        sig.add_hash(50);
        assert_eq!(sig.hashes, vec![100, 50]); // Order depends on implementation detail ignored here

        // Try adding more than capacity (current simple impl ignores)
        sig.add_hash(200);
        assert_eq!(sig.hashes.len(), 2); // Should not exceed capacity
    }

    #[test]
    fn test_jaccard_identical() {
        let sig1 = create_test_sig("sig1", vec![1, 2, 3, 4, 5]);
        let jaccard = sig1.jaccard(&sig1).unwrap();
        assert_relative_eq!(jaccard, 1.0);
    }

    #[test]
    fn test_jaccard_disjoint() {
        let sig1 = create_test_sig("sig1", vec![1, 2, 3]);
        let sig2 = create_test_sig("sig2", vec![4, 5, 6]);
        let jaccard = sig1.jaccard(&sig2).unwrap();
        assert_relative_eq!(jaccard, 0.0);
    }

    #[test]
    fn test_jaccard_partial_overlap() {
        let sig1 = create_test_sig("sig1", vec![1, 2, 3, 4]); // Set {1, 2, 3, 4}
        let sig2 = create_test_sig("sig2", vec![3, 4, 5, 6]); // Set {3, 4, 5, 6}
                                                              // Intersection {3, 4}, size = 2
                                                              // Union {1, 2, 3, 4, 5, 6}, size = 6
                                                              // Jaccard = 2 / 6 = 1/3
        let jaccard = sig1.jaccard(&sig2).unwrap();
        assert_relative_eq!(jaccard, 1.0 / 3.0, epsilon = 1e-9);
    }

    #[test]
    fn test_jaccard_subset() {
        let sig1 = create_test_sig("sig1", vec![1, 2, 3, 4]); // Set {1, 2, 3, 4}
        let sig2 = create_test_sig("sig2", vec![1, 2]); // Set {1, 2}
                                                        // Intersection {1, 2}, size = 2
                                                        // Union {1, 2, 3, 4}, size = 4
                                                        // Jaccard = 2 / 4 = 0.5
        let jaccard = sig1.jaccard(&sig2).unwrap();
        // Note: Jaccard might be less accurate if num_hashes differ significantly in real MinHash
        assert_relative_eq!(jaccard, 0.5, epsilon = 1e-9);
    }

    #[test]
    fn test_jaccard_empty() {
        let sig1 = create_test_sig("sig1", vec![1, 2]);
        let sig_empty = create_test_sig("empty", vec![]);
        assert_relative_eq!(sig1.jaccard(&sig_empty).unwrap(), 0.0);
        assert_relative_eq!(sig_empty.jaccard(&sig1).unwrap(), 0.0);
        assert_relative_eq!(sig_empty.jaccard(&sig_empty).unwrap(), 1.0); // Empty vs Empty
    }

    #[test]
    fn test_jaccard_incompatible() {
        let sig1 = create_test_sig("sig1", vec![1, 2]);
        let mut sig2 = create_test_sig("sig2", vec![3, 4]);
        sig2.kmer_size = 31; // Different k
        assert!(sig1.jaccard(&sig2).is_err());

        let mut sig3 = create_test_sig("sig3", vec![5, 6]);
        sig3.algorithm = "other_sketch".to_string(); // Different algorithm
        assert!(sig1.jaccard(&sig3).is_err());
    }

    // TODO: Add tests for merge() once implemented.
}
