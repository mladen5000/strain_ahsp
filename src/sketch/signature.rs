//! Genomic signature representations.
//!
//! This module provides implementations for different types of genomic signatures.

use bincode::{Decode, Encode};
use nthash::NtHashIterator;
use serde::{Deserialize, Serialize};
use std::collections::hash_map::DefaultHasher;
use std::collections::BinaryHeap; // Added for efficient intersection
use std::collections::HashSet; // Added for efficient intersection
use std::hash::{Hash, Hasher};
use std::path::{Path, PathBuf}; // Added Path for function args

// --- Generic Signature (Sketch) ---

/// Represents the core sketch data, typically a collection of hash values.
/// This struct is agnostic to the source of the items that were hashed (e.g., k-mers, words).
/// Its primary purpose is to hold the sketch and enable similarity comparisons
/// between sketches generated with compatible parameters.
#[derive(Debug, Clone, Serialize, Deserialize, Decode, Encode, PartialEq, Eq)]
pub struct Signature {
    // Algorithm used (e.g., "minhash", "scaled_minhash"). Crucial for interpreting hashes
    // and choosing the correct similarity calculation.
    pub algorithm: String,

    // The actual hash values comprising the sketch.
    // For MinHash, this is the list of minimum hashes.
    // For Scaled MinHash, this is hashes below a threshold.
    pub hashes: Vec<u64>,

    // --- Parameters defining the sketch resolution/size ---

    // Target number of hashes for fixed-size sketches (e.g., standard MinHash).
    // Set to 0 if using a variable-size method like scaled MinHash.
    pub num_hashes: usize,

    // Scaling factor for scaled MinHash (e.g., sourmash).
    // If > 0, the sketch contains hashes H where H < max_hash / scaled.
    // Set to 0 if not using scaled MinHash.
    pub scaled: u64,
    // Maximum possible hash value used by the hashing function. Needed for some
    // calculations, especially with scaled MinHash, though often implicit (e.g., u64::MAX).
    // Can be omitted if always using u64::MAX or if handled elsewhere.
    // pub max_hash: u64,
}

impl Signature {
    /// Creates a new, empty sketch container with given parameters.
    pub fn new(algorithm: String, num_hashes: usize, scaled: u64) -> Self {
        // Basic validation: Can't have both num_hashes > 0 and scaled > 0 easily defined
        // This depends on the specific algorithm implementation details.
        // assert!(!(num_hashes > 0 && scaled > 0), "Cannot specify both num_hashes and scaled factor for a single sketch");

        Signature {
            algorithm,
            // Pre-allocate if using fixed-size num_hashes for efficiency
            hashes: if num_hashes > 0 {
                Vec::with_capacity(num_hashes)
            } else {
                Vec::new()
            },
            num_hashes,
            scaled,
            // max_hash: u64::MAX, // Example default
        }
    }

    /// Checks if the sketch is empty.
    pub fn is_empty(&self) -> bool {
        self.hashes.is_empty()
    }

    /// Returns the effective size of the sketch.
    pub fn size(&self) -> usize {
        self.hashes.len()
    }

    /// Calculates the Jaccard similarity estimate between this sketch and another.
    /// Assumes a MinHash-like sketch (standard or scaled).
    ///
    /// # Arguments
    /// * `other` - Another Signature sketch to compare with.
    ///
    /// # Returns
    /// The Jaccard similarity estimate (between 0.0 and 1.0), or None if
    /// sketches are incompatible (different algorithms, incompatible parameters).
    pub fn estimate_jaccard(&self, other: &Signature) -> Option<f64> {
        // --- Compatibility Checks ---
        if self.algorithm != other.algorithm {
            return None; // Different algorithms cannot be compared directly
        }

        // Check if parameters make sense for comparison
        // For standard MinHash (num_hashes > 0), num_hashes should ideally match,
        // but comparison using the minimum is common.
        // For scaled MinHash (scaled > 0), scaled factors must match.
        if self.scaled > 0 {
            if self.scaled != other.scaled {
                return None; // Scaled MinHash requires matching scaling factors
            }
        } else if self.num_hashes > 0 {
            // For fixed num_hashes, we can compare even if nums differ, using the smaller num.
            // Warning: This estimate is less reliable if num_hashes differ significantly.
            if other.num_hashes == 0 {
                return None;
            } // Cannot compare fixed num with scaled (scaled=0)
        } else {
            // Both num_hashes and scaled are 0 - undefined sketch type? Or empty?
            if other.num_hashes > 0 || other.scaled > 0 {
                return None;
            } // Incompatible
              // If both are 0 and empty, similarity is arguably 1.0? Or 0.0? Or undefined? Let's return 0.0 for empty.
        }

        if self.is_empty() || other.is_empty() {
            // Jaccard of empty sets is often considered 0, or 1 if both are empty.
            // Let's return 0.0 if either is empty (unless both are empty AND compatible)
            if self.is_empty()
                && other.is_empty()
                && self.algorithm == other.algorithm
                && self.scaled == other.scaled
                && self.num_hashes == other.num_hashes
            {
                return Some(1.0); // J(empty, empty) = 1
            } else {
                return Some(0.0); // J(A, empty) = 0 if A is not empty
            }
        }

        // --- Calculate Intersection ---
        // Use HashSet for efficiency with larger sketches
        let self_hashes: HashSet<u64> = self.hashes.iter().cloned().collect();
        let mut intersection_size = 0;
        for hash in &other.hashes {
            if self_hashes.contains(hash) {
                intersection_size += 1;
            }
        }

        // --- Estimate Jaccard based on algorithm type ---
        if self.scaled > 0 {
            // Scaled MinHash: J = |Intersection| / |Union|
            // Estimate Union size = total unique hashes observed across both sketches
            // This simple approximation uses len(A) + len(B) - intersection, which is exact if hashes are unique within each vec
            let union_size = self.hashes.len() + other.hashes.len() - intersection_size;
            if union_size == 0 {
                return Some(1.0); // Should be caught by is_empty checks, but safe fallback
            }
            Some(intersection_size as f64 / union_size as f64)
        } else if self.num_hashes > 0 {
            // Standard MinHash: J approx |Intersection| / num_hashes (or min(num_hashes))
            let min_num_hashes = self.num_hashes.min(other.num_hashes);
            if min_num_hashes == 0 {
                return None;
            } // Should not happen if num_hashes > 0 check passed
            Some(intersection_size as f64 / min_num_hashes as f64)
        } else {
            // Undefined case (both num_hashes and scaled are 0)
            None // Or handle as appropriate if this state is valid
        }
    }
}

impl Default for Signature {
    fn default() -> Self {
        Self {
            algorithm: "minhash".to_string(),
            hashes: Vec::new(),
            num_hashes: 0,
            scaled: 0,
            // max_hash: u64::MAX,
        }
    }
}

// --- Kmer Signature (Sequence Sketch Context) ---

/// Represents a signature/sketch derived specifically from the k-mers of a sequence.
/// It bundles the actual sketch (`Signature`) with essential contextual information
/// like k-mer size and molecule type required for meaningful interpretation and comparison.
#[derive(Debug, Clone, Serialize, Deserialize, Decode, Encode)]
pub struct KmerSignature {
    // The core sketch data generated from the k-mers
    pub sketch: Signature,

    // K-mer size used to generate the sketch
    pub kmer_size: usize,

    // Type of molecule the sequence represents (e.g., "DNA", "protein")
    pub molecule_type: String,

    // Optional name for the signature (e.g., sequence ID)
    pub name: Option<String>,

    // Optional filename and path information
    pub filename: Option<String>,
    pub path: Option<PathBuf>,
}

impl KmerSignature {
    pub fn is_initialized(&self) -> bool {
        // Check relevant fields that indicate proper initialization
        self.kmer_size > 0 && self.molecule_type.len() > 0
    }

    /// Calculates the Jaccard similarity between this KmerSignature and another.
    /// Ensures that k-mer sizes and molecule types are compatible before comparing sketches.
    ///
    /// # Arguments
    /// * `other` - Another KmerSignature to compare with.
    ///
    /// # Returns
    /// The Jaccard similarity estimate (0.0 to 1.0) if comparable, otherwise None.
    pub fn jaccard_similarity(&self, other: &KmerSignature) -> Option<f64> {
        // Check context first
        if self.kmer_size != other.kmer_size {
            return None;
        }
        // Molecule types should be compatible for meaningful comparison
        if !self.are_molecule_types_compatible(&other.molecule_type) {
            return None;
        }

        // Delegate to the underlying Signature's jaccard_similarity
        Some(self.sketch.estimate_jaccard(&other.sketch).unwrap_or(0.0))
    }

    /// Checks if molecule types are compatible for comparison
    fn are_molecule_types_compatible(&self, other_type: &str) -> bool {
        // DNA and RNA can be compared (they use same canonical k-mers)
        let is_dna = self.molecule_type.eq_ignore_ascii_case("DNA");
        let is_rna = self.molecule_type.eq_ignore_ascii_case("RNA");
        let other_is_dna = other_type.eq_ignore_ascii_case("DNA");
        let other_is_rna = other_type.eq_ignore_ascii_case("RNA");

        (is_dna || is_rna) && (other_is_dna || other_is_rna)
            || self.molecule_type.eq_ignore_ascii_case(other_type)
    }

    /// Adds a sequence to the signature by processing its k-mers and updating the sketch.
    /// Uses ntHash for hashing. If molecule_type is "DNA" or "RNA" (case-insensitive),
    /// it processes canonical k-mer hashes.
    ///
    /// Returns an error if the sequence is invalid, k-mer size is incompatible,
    /// or hashing/sketching fails.
    pub fn add_sequence(&mut self, sequence: &[u8]) -> Result<(), String> {
        // Determine if we should use canonical k-mers
        let use_canonical = self.molecule_type.eq_ignore_ascii_case("DNA")
            || self.molecule_type.eq_ignore_ascii_case("RNA");

        // Create ntHash Iterator
        let hasher = NtHashIterator::new(sequence, self.kmer_size)
            .map_err(|_| format!("ntHash failed to initialize for k={}", self.kmer_size))?;

        // Process k-mer hashes based on sketch type (fixed-size MinHash vs scaled MinHash)
        if self.sketch.num_hashes > 0 {
            // Fixed-size MinHash: Keep the smallest num_hashes unique values
            let mut heap = BinaryHeap::from(self.sketch.hashes.clone());

            for hash_value in hasher {
                let canonical_hash = if use_canonical {
                    // For DNA/RNA, hash both k-mer and its reverse complement, take the smaller value
                    let rc_hash = hash_value.rotate_left(1); // Simple way to get a different hash for RC
                    hash_value.min(rc_hash)
                } else {
                    hash_value
                };

                if heap.len() < self.sketch.num_hashes {
                    heap.push(canonical_hash);
                } else if let Some(&max_hash) = heap.peek() {
                    if canonical_hash < max_hash {
                        heap.pop();
                        heap.push(canonical_hash);
                    }
                }
            }

            self.sketch.hashes = heap.into_sorted_vec();
        } else if self.sketch.scaled > 0 {
            // Scaled MinHash: Keep all hashes below the threshold
            let threshold = u64::MAX / self.sketch.scaled;
            let mut kept_hashes = HashSet::new();

            for hash_value in hasher {
                let canonical_hash = if use_canonical {
                    let rc_hash = hash_value.rotate_left(1);
                    hash_value.min(rc_hash)
                } else {
                    hash_value
                };

                if canonical_hash < threshold {
                    kept_hashes.insert(canonical_hash);
                }
            }

            self.sketch.hashes = kept_hashes.into_iter().collect();
            self.sketch.hashes.sort_unstable();
        } else {
            return Err(format!(
                "Invalid sketch parameters: num_hashes={}, scaled={}",
                self.sketch.num_hashes, self.sketch.scaled
            ));
        }

        Ok(())
    }
}

// --- Multi Resolution Signature ---

/// Resolution level for hierarchical sketches (Conceptual).
#[derive(Debug, Clone, Hash, PartialEq, Eq, Serialize, Deserialize)]
pub enum ResolutionLevel {
    Macro,      // Coarse resolution (e.g., smaller k, larger scale/fewer hashes)
    Meso,       // Medium resolution
    Micro,      // Fine resolution (e.g., larger k, smaller scale/more hashes)
    Custom(u8), // Custom resolution identifier
}

/// A multi-resolution genomic signature, holding several KmerSignatures
/// likely generated with different parameters (k-mer size, sketch size)
/// to capture similarity at different scales.
#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, Default)]
pub struct MultiResolutionSignature {
    // Identifier for the taxon/genome this signature represents.
    pub taxon_id: String,
    // Taxonomic lineage information.
    pub lineage: Vec<String>,
    // The collection of signatures at different resolutions.
    // Use a HashMap or Vec with associated ResolutionLevel if needed for lookup.
    // Using Vec implies order matters (e.g., index 0=Macro, 1=Meso, 2=Micro).
    pub levels: Vec<KmerSignature>,
    // We remove the redundant fields 'macro_signature', etc.
    // They should be accessed via the `levels` vector.
    // #[serde(skip)] pub macro_signature: KmerSignature, // Removed
    // #[serde(skip)] pub meso_signature: KmerSignature,  // Removed
    // #[serde(skip)] pub micro_signature: KmerSignature, // Removed
}

impl MultiResolutionSignature {
    /// Create a new multi-resolution signature container.
    pub fn new(taxon_id: String, lineage: Vec<String>) -> Self {
        MultiResolutionSignature {
            taxon_id,
            lineage,
            levels: Vec::new(), // Initialize levels vector
        }
    }

    /// Adds a KmerSignature for a specific resolution level.
    /// Note: This simple version just adds to the Vec. A real implementation
    /// might associate it with a ResolutionLevel enum or ensure specific ordering.
    pub fn add_level(&mut self, signature: KmerSignature) {
        self.levels.push(signature);
    }

    /// Calculate similarity between this signature and another
    pub fn similarity(&self, other: &Self, weights: Option<Vec<f64>>) -> Option<f64> {
        if self.levels.is_empty() || other.levels.is_empty() {
            return None;
        }

        // Use equal weights if none provided
        let num_levels = self.levels.len().min(other.levels.len());
        let weights = weights.unwrap_or_else(|| {
            let weight = 1.0 / num_levels as f64;
            vec![weight; num_levels]
        });

        let mut total_similarity = 0.0;
        for (i, (self_level, other_level)) in self
            .levels
            .iter()
            .zip(other.levels.iter())
            .take(num_levels)
            .enumerate()
        {
            if let Some(sim) = self_level.jaccard_similarity(other_level) {
                total_similarity += weights[i] * sim;
            } else {
                return None; // Unable to compare signatures at this level
            }
        }

        Some(total_similarity)
    }
}

// --- Builder Pattern ---

/// Builder for creating KmerSignature objects more fluently.
#[derive(Debug, Clone)]
pub struct KmerSignatureBuilder {
    kmer_size: usize,
    molecule_type: String,
    algorithm: String,
    num_hashes: usize,
    scaled: u64,
    name: Option<String>,
    path: Option<PathBuf>,
}

impl KmerSignatureBuilder {
    /// Creates a new builder with required parameters.
    pub fn new(
        kmer_size: usize,
        molecule_type: &str,
        algorithm: &str,
        num_hashes: usize,
        scaled: u64,
    ) -> Self {
        KmerSignatureBuilder {
            kmer_size,
            molecule_type: molecule_type.to_string(),
            algorithm: algorithm.to_string(),
            num_hashes,
            scaled,
            name: None,
            path: None,
        }
    }

    /// Sets the optional name for the signature.
    pub fn name(mut self, name: &str) -> Self {
        self.name = Some(name.to_string());
        self
    }

    /// Sets the optional source path for the signature.
    pub fn source(mut self, path: &Path) -> Self {
        self.path = Some(path.to_path_buf());
        self
    }

    /// Builds the KmerSignature.
    pub fn build(&self) -> KmerSignature {
        let mut signature = KmerSignature {
            sketch: Signature::new(self.algorithm.clone(), self.num_hashes, self.scaled),
            kmer_size: self.kmer_size,
            molecule_type: self.molecule_type.clone(),
            name: self.name.clone(),
            filename: self
                .path
                .as_ref()
                .and_then(|p| p.file_name())
                .map(|n| n.to_string_lossy().into_owned()),
            path: self.path.clone(),
        };
        signature
    }
}

// --- Deprecated/Removed ---
// The old KmerSignature struct and Sketch trait are removed as their
// functionality is either integrated into the new KmerSignature or wasn't
// representing a sketch in the MinHash sense.

// --- Tests ---
#[cfg(test)]
mod tests {
    use super::*;

    // Helper to create a basic KmerSignature for testing
    fn create_test_kmer_sig(name: &str, k: usize, num: usize, hashes: Vec<u64>) -> KmerSignature {
        let mut ksig = KmerSignatureBuilder::new(k, "DNA", "minhash", num, 0)
            .name(name)
            .build();
        ksig.sketch.hashes = hashes;
        ksig
    }
    // Helper to create a scaled KmerSignature for testing
    fn create_scaled_test_kmer_sig(
        name: &str,
        k: usize,
        scale: u64,
        hashes: Vec<u64>,
    ) -> KmerSignature {
        let mut ksig = KmerSignatureBuilder::new(k, "DNA", "scaled_minhash", 0, scale)
            .name(name)
            .build();
        ksig.sketch.hashes = hashes;
        ksig
    }

    #[test]
    fn test_jaccard_minhash_equal() {
        let sig1 = create_test_kmer_sig("sig1", 21, 5, vec![1, 2, 3, 4, 5]);
        let sig2 = create_test_kmer_sig("sig2", 21, 5, vec![1, 2, 3, 4, 5]);
        assert_eq!(sig1.sketch.estimate_jaccard(&sig2.sketch), Some(1.0));
        assert_eq!(sig2.sketch.estimate_jaccard(&sig1.sketch), Some(1.0)); // Symmetric
    }

    #[test]
    fn test_jaccard_scaled_equal() {
        let sig1 = create_scaled_test_kmer_sig("sig1", 21, 1000, vec![10, 20, 30]);
        let sig2 = create_scaled_test_kmer_sig("sig2", 21, 1000, vec![10, 20, 30]);
        // Intersection=3, Union=3+3-3=3. J=3/3=1.0
        assert_eq!(sig1.sketch.estimate_jaccard(&sig2.sketch), Some(1.0));
    }

    #[test]
    fn test_jaccard_minhash_different() {
        let sig1 = create_test_kmer_sig("sig1", 21, 5, vec![1, 2, 3, 4, 5]);
        let sig2 = create_test_kmer_sig("sig2", 21, 5, vec![6, 7, 8, 9, 10]);
        // Intersection=0. J=0/5=0.0
        assert_eq!(sig1.sketch.estimate_jaccard(&sig2.sketch), Some(0.0));
    }

    #[test]
    fn test_jaccard_scaled_different() {
        let sig1 = create_scaled_test_kmer_sig("sig1", 21, 1000, vec![10, 20, 30]);
        let sig2 = create_scaled_test_kmer_sig("sig2", 21, 1000, vec![40, 50, 60]);
        // Intersection=0, Union=3+3-0=6. J=0/6=0.0
        assert_eq!(sig1.sketch.estimate_jaccard(&sig2.sketch), Some(0.0));
    }

    #[test]
    fn test_jaccard_minhash_partial() {
        let sig1 = create_test_kmer_sig("sig1", 21, 5, vec![1, 2, 3, 4, 5]);
        let sig2 = create_test_kmer_sig("sig2", 21, 5, vec![1, 2, 3, 9, 10]);
        // Intersection=3. J=3/5=0.6
        assert_eq!(sig1.sketch.estimate_jaccard(&sig2.sketch), Some(0.6));
    }

    #[test]
    fn test_jaccard_scaled_partial() {
        let sig1 = create_scaled_test_kmer_sig("sig1", 21, 1000, vec![10, 20, 30, 40, 50]); // size 5
        let sig2 = create_scaled_test_kmer_sig("sig2", 21, 1000, vec![10, 20, 30, 60, 70]); // size 5
                                                                                            // Intersection=3, Union=5+5-3=7. J=3/7
        assert!((sig1.sketch.estimate_jaccard(&sig2.sketch).unwrap() - (3.0 / 7.0)).abs() < 1e-9);
    }

    #[test]
    fn test_jaccard_minhash_different_num_hashes() {
        let sig1 = create_test_kmer_sig("sig1", 21, 10, vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
        let sig2 = create_test_kmer_sig("sig2", 21, 5, vec![1, 2, 3, 11, 12]);
        // Intersection=3. min_num_hashes=5. J=3/5=0.6
        assert_eq!(sig1.sketch.estimate_jaccard(&sig2.sketch), Some(0.6));
    }

    #[test]
    fn test_jaccard_incompatible_k() {
        let sig1 = create_test_kmer_sig("sig1", 21, 5, vec![1, 2, 3, 4, 5]);
        let sig2 = create_test_kmer_sig("sig2", 31, 5, vec![1, 2, 3, 4, 5]); // Different k
        assert_eq!(sig1.sketch.estimate_jaccard(&sig2.sketch), None);
    }

    #[test]
    fn test_jaccard_incompatible_scaled() {
        let sig1 = create_scaled_test_kmer_sig("sig1", 21, 1000, vec![10, 20, 30]);
        let sig2 = create_scaled_test_kmer_sig("sig2", 21, 2000, vec![10, 20, 30]); // Different scale
        assert_eq!(sig1.sketch.estimate_jaccard(&sig2.sketch), None);
    }

    #[test]
    fn test_jaccard_incompatible_algo() {
        let sig1 = create_test_kmer_sig("sig1", 21, 5, vec![1, 2, 3, 4, 5]);
        let mut sig2 = create_test_kmer_sig("sig2", 21, 5, vec![1, 2, 3, 4, 5]);
        sig2.sketch.algorithm = "other_algo".to_string();
        assert_eq!(sig1.sketch.estimate_jaccard(&sig2.sketch), None);
    }

    #[test]
    fn test_kmer_signature_builder() {
        let builder = KmerSignatureBuilder::new(21, "DNA", "minhash", 500, 0)
            .name("MyGenome")
            .source(Path::new("/data/mygenome.fna"));

        let ksig = builder.build();

        assert_eq!(ksig.kmer_size, 21);
        assert_eq!(ksig.molecule_type, "DNA");
        assert_eq!(ksig.sketch.algorithm, "minhash");
        assert_eq!(ksig.sketch.num_hashes, 500);
        assert_eq!(ksig.sketch.scaled, 0);
        assert_eq!(ksig.name.as_deref(), Some("MyGenome"));
        assert_eq!(ksig.path.as_deref(), Some(Path::new("/data/mygenome.fna")));
        assert_eq!(ksig.filename.as_deref(), Some("mygenome.fna"));
        assert!(ksig.sketch.hashes.is_empty()); // Builder doesn't add sequences
    }

    #[test]
    fn test_multi_resolution_similarity() {
        let mut mrs1 = MultiResolutionSignature::new("tax1".to_string(), vec![]);
        let mut mrs2 = MultiResolutionSignature::new("tax2".to_string(), vec![]);

        // Create test signatures for different resolution levels
        // Level 1 (Macro) - 5/10 overlap
        let mut level1_sig1 = KmerSignatureBuilder::new(15, "DNA", "minhash", 10, 0)
            .name("l1s1")
            .build();
        level1_sig1.sketch.hashes = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];

        let mut level1_sig2 = KmerSignatureBuilder::new(15, "DNA", "minhash", 10, 0)
            .name("l1s2")
            .build();
        level1_sig2.sketch.hashes = vec![1, 2, 3, 4, 5, 11, 12, 13, 14, 15];

        // Level 2 (Meso) - 2/20 overlap (using scaled MinHash)
        let mut level2_sig1 = KmerSignatureBuilder::new(11, "DNA", "scaled_minhash", 0, 20)
            .name("l2s1")
            .build();
        level2_sig1.sketch.hashes = vec![10, 20, 30, 40, 50];

        let mut level2_sig2 = KmerSignatureBuilder::new(11, "DNA", "scaled_minhash", 0, 20)
            .name("l2s2")
            .build();
        level2_sig2.sketch.hashes = vec![10, 20, 60, 70, 80];

        // Add levels to signatures
        mrs1.add_level(level1_sig1);
        mrs1.add_level(level2_sig1);
        mrs2.add_level(level1_sig2);
        mrs2.add_level(level2_sig2);

        // Test with default weights (equal weighting)
        let sim_default = mrs1.similarity(&mrs2, None);
        assert!(sim_default.is_some());
        // Expected similarity: (0.5 + 0.4) / 2 = 0.45
        // Level 1: 5 shared / 10 total = 0.5
        // Level 2: 2 shared / 5 total = 0.4
        assert!((sim_default.unwrap() - 0.45).abs() < 1e-9);

        // Test with custom weights
        let weights = vec![0.3, 0.7];
        let sim_custom = mrs1.similarity(&mrs2, Some(weights));
        assert!(sim_custom.is_some());
        // Expected: (0.3 * 0.5) + (0.7 * 0.4) = 0.15 + 0.28 = 0.43
        assert!((sim_custom.unwrap() - 0.43).abs() < 1e-9);
    }
}
