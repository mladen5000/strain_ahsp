//! Genomic signature representations.
//!
//! This module provides implementations for different types of genomic signatures.

use bincode::{Decode, Encode};
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

/// Represents a sketch/signature of a sequence, typically a collection of hash values.
/// For example, a MinHash signature would be a set of smallest hash values.
#[derive(Debug, Clone, Serialize, Deserialize, Decode, Encode)]
pub struct Signature {
    pub algorithm: String, // Algorithm used to generate the signature (e.g., "minhash")
    pub kmer_size: usize,  // K-mer size used
    pub num_hashes: usize, // Number of hashes requested (actual might be less)
    pub name: Option<String>, // Optional name for the signature (e.g., sequence ID)
    pub filename: Option<String>, // Filename source, if sketched from a file
    pub path: Option<PathBuf>, // File path, if sketched from a file
    pub hashes: Vec<u64>,  // The actual hash values comprising the signature
}

impl Signature {
    /// Creates a new signature with the given parameters.
    pub fn new(algorithm: String, kmer_size: usize, num_hashes: usize) -> Self {
        Signature {
            algorithm,
            kmer_size,
            num_hashes,
            name: None,
            filename: None,
            path: None,
            hashes: Vec::with_capacity(num_hashes), // Pre-allocate for efficiency
        }
    }

    /// Creates a new signature with a specific name.
    pub fn with_name(algorithm: String, kmer_size: usize, num_hashes: usize, name: String) -> Self {
        Signature {
            algorithm,
            kmer_size,
            num_hashes,
            name: Some(name),
            filename: None,
            path: None,
            hashes: Vec::with_capacity(num_hashes),
        }
    }

    /// Calculates the Jaccard similarity between this signature and another.
    ///
    /// # Arguments
    ///
    /// * `other` - Another Signature to compare with.
    ///
    /// # Returns
    ///
    /// The Jaccard similarity estimate (between 0.0 and 1.0).
    pub fn jaccard_similarity(&self, other: &Signature) -> f64 {
        if self.algorithm != other.algorithm || self.kmer_size != other.kmer_size {
            // If parameters don't match, similarity is meaningless
            // Could return 0.0, or error - depends on use case
            return 0.0;
        }

        // For signatures with different num_hashes, take the smaller one
        let min_num_hashes = self.num_hashes.min(other.num_hashes);

        if min_num_hashes == 0 || self.hashes.is_empty() || other.hashes.is_empty() {
            return 0.0; // No valid comparison possible
        }

        // Calculate intersection size (for MinHash & similar techniques)
        // This is a simple set intersection approach. More sophisticated algorithms might
        // compute the intersection differently depending on the technique.
        let mut intersection_size = 0;

        // Simple O(n^2) approach for small signatures
        if self.hashes.len() <= 100 && other.hashes.len() <= 100 {
            for hash in &self.hashes {
                if other.hashes.contains(hash) {
                    intersection_size += 1;
                }
            }
        } else {
            // For larger signatures, use a set-based approach
            let self_set: std::collections::HashSet<&u64> = self.hashes.iter().collect();
            for hash in &other.hashes {
                if self_set.contains(hash) {
                    intersection_size += 1;
                }
            }
        }

        // For standard MinHash, Jaccard similarity is estimated as:
        // |intersection| / min_num_hashes
        intersection_size as f64 / min_num_hashes as f64
    }
}

impl Default for Signature {
    fn default() -> Self {
        Self {
            algorithm: "empty".to_string(),
            kmer_size: 0,
            num_hashes: 0,
            name: None,
            filename: None,
            path: None,
            hashes: Vec::new(),
        }
    }
}

/// Resolution level for hierarchical sketches.
#[derive(Debug, Clone, Hash, PartialEq, Eq, Serialize, Deserialize)]
pub enum ResolutionLevel {
    Macro,      // Coarse resolution
    Meso,       // Medium resolution
    Micro,      // Fine resolution
    Custom(u8), // Custom resolution with kmer size
}

/// A multi-resolution genomic signature.
#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, Default)]
pub struct MultiResolutionSignature {
    pub taxon_id: String,
    pub lineage: Vec<String>,
    pub levels: Vec<Signature>,

    // These are added for compatibility with existing code
    #[serde(skip)]
    pub macro_signature: Signature,
    #[serde(skip)]
    pub meso_signature: Signature,
    #[serde(skip)]
    pub micro_signature: Signature,
}

impl MultiResolutionSignature {
    /// Create a new multi-resolution signature
    pub fn new(taxon_id: String, lineage: Vec<String>) -> Self {
        MultiResolutionSignature {
            taxon_id,
            lineage,
            levels: Vec::new(),
            macro_signature: Signature::default(),
            meso_signature: Signature::default(),
            micro_signature: Signature::default(),
        }
    }

    /// Calculate similarity between this signature and another
    pub fn similarity(&self, other: &Self, weights: Option<Vec<f64>>) -> f64 {
        let default_weights = vec![0.2, 0.3, 0.5]; // Example weights for macro, meso, micro
        let weights = weights.unwrap_or(default_weights);

        // Example weighted combination of similarities
        let macro_sim = self
            .macro_signature
            .jaccard_similarity(&other.macro_signature);
        let meso_sim = self
            .meso_signature
            .jaccard_similarity(&other.meso_signature);
        let micro_sim = self
            .micro_signature
            .jaccard_similarity(&other.micro_signature);

        weights[0] * macro_sim + weights[1] * meso_sim + weights[2] * micro_sim
    }
}

/// Represents a k-mer signature, used for sequences.
#[derive(Debug, Clone, Serialize, Deserialize, Decode, Encode)]
pub struct KmerSignature {
    pub kmers: Vec<Vec<u8>>,
    pub counts: Vec<u32>,
    pub total_kmers: u64,
    pub(crate) sketch_size: usize,
    pub sketch: Vec<u64>,
}

impl KmerSignature {
    pub fn is_initialized(&self) -> bool {
        // Check relevant fields that indicate proper initialization
        // For example, check if sketch size > 0 or if there's data in the signature
        !self.sketch.is_empty() // Assuming there's a sketch field
    }
}

/// Builder for creating signature objects.
#[derive(Debug)]
pub struct SignatureBuilder {
    pub algorithm: String,
    pub kmer_size: usize,
    pub num_hashes: usize,
}

impl SignatureBuilder {
    /// Creates a new signature builder.
    pub fn new(algorithm: &str, kmer_size: usize, num_hashes: usize) -> Self {
        SignatureBuilder {
            algorithm: algorithm.to_string(),
            kmer_size,
            num_hashes,
        }
    }

    /// Builds a new empty signature.
    pub fn build(&self) -> Signature {
        Signature::new(self.algorithm.clone(), self.kmer_size, self.num_hashes)
    }

    /// Builds a new signature with the given name.
    pub fn build_with_name(&self, name: &str) -> Signature {
        Signature::with_name(
            self.algorithm.clone(),
            self.kmer_size,
            self.num_hashes,
            name.to_string(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_jaccard_similarity_equal() {
        let mut sig1 = Signature::new("minhash".to_string(), 21, 5);
        let mut sig2 = Signature::new("minhash".to_string(), 21, 5);

        // Identical hashes
        sig1.hashes = vec![1, 2, 3, 4, 5];
        sig2.hashes = vec![1, 2, 3, 4, 5];

        assert_eq!(sig1.jaccard_similarity(&sig2), 1.0);
        assert_eq!(sig2.jaccard_similarity(&sig1), 1.0); // Should be symmetric
    }

    #[test]
    fn test_jaccard_similarity_different() {
        let mut sig1 = Signature::new("minhash".to_string(), 21, 5);
        let mut sig2 = Signature::new("minhash".to_string(), 21, 5);

        // Completely different hashes
        sig1.hashes = vec![1, 2, 3, 4, 5];
        sig2.hashes = vec![6, 7, 8, 9, 10];

        assert_eq!(sig1.jaccard_similarity(&sig2), 0.0);
    }

    #[test]
    fn test_jaccard_similarity_partial() {
        let mut sig1 = Signature::new("minhash".to_string(), 21, 5);
        let mut sig2 = Signature::new("minhash".to_string(), 21, 5);

        // 3/5 shared hashes
        sig1.hashes = vec![1, 2, 3, 4, 5];
        sig2.hashes = vec![1, 2, 3, 9, 10];

        assert_eq!(sig1.jaccard_similarity(&sig2), 0.6);
    }

    #[test]
    fn test_jaccard_similarity_different_sizes() {
        let mut sig1 = Signature::new("minhash".to_string(), 21, 10);
        let mut sig2 = Signature::new("minhash".to_string(), 21, 5);

        // First has 10 hashes, second has 5, with 3 shared
        sig1.hashes = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        sig2.hashes = vec![1, 2, 3, 11, 12];

        // Since we're using min_num_hashes as denominator, this is 3/5
        assert_eq!(sig1.jaccard_similarity(&sig2), 0.6);
    }

    #[test]
    fn test_signature_builder() {
        let builder = SignatureBuilder::new("minhash", 21, 5);
        let sig = builder.build();

        assert_eq!(sig.algorithm, "minhash");
        assert_eq!(sig.kmer_size, 21);
        assert_eq!(sig.num_hashes, 5);
        assert!(sig.hashes.is_empty());
    }
}
