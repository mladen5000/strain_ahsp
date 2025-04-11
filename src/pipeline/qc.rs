use needletail::{parse_fastx_file, Sequence};
use nthash::NtHashIterator;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::path::Path;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum SignatureError {
    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),
    #[error("Invalid sequence format: {0}")]
    SequenceFormatError(String),
    #[error("Invalid k-mer size: {0}")]
    InvalidKmerSize(usize),
}

/// Resolution levels for genomic signatures
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ResolutionLevel {
    Macro, // Species level (larger k-mers, e.g., 31)
    Meso,  // Strain group level (mid-size k-mers, e.g., 21)
    Micro, // Individual strain level (SNVs and smaller k-mers, e.g., 11)
}

/// K-mer signature using MinHash sketch
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KmerSignature {
    /// K-mer size
    pub k: usize,
    /// Sketch size
    pub sketch_size: usize,
    /// MinHash sketch values
    pub sketch: Vec<u64>,
    /// Total k-mers processed (for abundance estimation)
    pub total_kmers: usize,
}

impl KmerSignature {
    /// Create a new empty k-mer signature
    pub fn new(k: usize, sketch_size: usize) -> Result<Self, SignatureError> {
        if k < 3 || k > 31 {
            return Err(SignatureError::InvalidKmerSize(k));
        }

        Ok(KmerSignature {
            k,
            sketch_size,
            sketch: Vec::with_capacity(sketch_size),
            total_kmers: 0,
        })
    }

    /// Add sequence to the signature
    pub fn add_sequence(&mut self, sequence: &[u8]) -> Result<(), SignatureError> {
        // Check sequence validity (ensuring ACGT only)
        if !sequence
            .iter()
            .all(|&b| b == b'A' || b == b'C' || b == b'G' || b == b'T' || b == b'N')
        {
            return Err(SignatureError::SequenceFormatError(
                "Sequence contains invalid characters".to_string(),
            ));
        }

        // Use ntHash to efficiently generate canonical k-mer hashes
        if let Ok(iterator) = NtHashIterator::new(sequence, (self.k as u8).into()) {
            // Update bottom sketch using MinHash technique
            for hash in iterator {
                self.total_kmers += 1;

                // Bottom-k MinHash implementation
                if self.sketch.len() < self.sketch_size {
                    // Sketch not full yet, add and sort
                    self.sketch.push(hash);
                    self.sketch.sort_unstable();
                } else if hash < self.sketch[self.sketch_size - 1] {
                    // Replace largest element and re-sort
                    self.sketch[self.sketch_size - 1] = hash;
                    self.sketch.sort_unstable();
                }
            }
        }

        Ok(())
    }

    /// Calculate Jaccard similarity between this signature and another
    pub fn jaccard_similarity(&self, other: &KmerSignature) -> f64 {
        if self.k != other.k {
            return 0.0; // Different k-mer sizes aren't comparable
        }

        let mut intersection_size = 0;
        let mut i = 0;
        let mut j = 0;

        // Count intersection using sorted sketches
        while i < self.sketch.len() && j < other.sketch.len() {
            if self.sketch[i] == other.sketch[j] {
                intersection_size += 1;
                i += 1;
                j += 1;
            } else if self.sketch[i] < other.sketch[j] {
                i += 1;
            } else {
                j += 1;
            }
        }

        // Calculate Jaccard index
        if self.sketch.is_empty() && other.sketch.is_empty() {
            return 1.0;
        }

        let union_estimate = self.sketch.len() + other.sketch.len() - intersection_size;
        intersection_size as f64 / union_estimate as f64
    }
}

/// Variant profile for fine-grained strain differentiation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VariantProfile {
    /// Maps position to observed nucleotide
    pub variants: HashMap<usize, u8>,
    /// Coverage at each position
    pub coverage: HashMap<usize, usize>,
}

impl VariantProfile {
    /// Create a new empty variant profile
    pub fn new() -> Self {
        VariantProfile {
            variants: HashMap::new(),
            coverage: HashMap::new(),
        }
    }

    /// Add a variant at a specific position
    pub fn add_variant(&mut self, position: usize, nucleotide: u8, count: usize) {
        // Update variant if new count is higher
        if let Some(current_coverage) = self.coverage.get(&position) {
            if count > *current_coverage {
                self.variants.insert(position, nucleotide);
                self.coverage.insert(position, count);
            }
        } else {
            self.variants.insert(position, nucleotide);
            self.coverage.insert(position, count);
        }
    }

    /// Calculate similarity between this variant profile and another
    pub fn similarity(&self, other: &VariantProfile) -> f64 {
        let mut shared_positions = 0;
        let mut matching_variants = 0;

        // Find positions covered in both profiles
        for position in self.variants.keys() {
            if other.variants.contains_key(position) {
                shared_positions += 1;

                if self.variants[position] == other.variants[position] {
                    matching_variants += 1;
                }
            }
        }

        if shared_positions == 0 {
            return 0.0;
        }

        matching_variants as f64 / shared_positions as f64
    }
}

/// Multi-resolution genomic signature
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MultiResolutionSignature {
    /// Taxonomic ID associated with this signature
    pub taxon_id: String,

    /// Full taxonomic lineage
    pub lineage: Vec<String>,

    /// Macro-level signature (species) - larger k-mers for uniqueness
    pub macro_signature: KmerSignature,

    /// Meso-level signature (strain group) - focused on variable regions
    pub meso_signature: KmerSignature,

    /// Micro-level signature (strain-specific variants)
    pub micro_signature: VariantProfile,

    /// Discriminatory power of each component
    pub weights: HashMap<ResolutionLevel, f64>,
}

impl MultiResolutionSignature {
    /// Create a new multi-resolution signature
    pub fn new(
        taxon_id: &str,
        lineage: Vec<String>,
        macro_k: usize,
        meso_k: usize,
        macro_sketch_size: usize,
        meso_sketch_size: usize,
    ) -> Result<Self, SignatureError> {
        Ok(MultiResolutionSignature {
            taxon_id: taxon_id.to_string(),
            lineage,
            macro_signature: KmerSignature::new(macro_k, macro_sketch_size)?,
            meso_signature: KmerSignature::new(meso_k, meso_sketch_size)?,
            micro_signature: VariantProfile::new(),
            weights: HashMap::new(),
        })
    }

    /// Add sequence to all resolution levels
    pub fn add_sequence(&mut self, sequence: &[u8]) -> Result<(), SignatureError> {
        // Update k-mer signatures
        self.macro_signature.add_sequence(sequence)?;
        self.meso_signature.add_sequence(sequence)?;

        // For micro-level (SNVs), we'd typically need alignment information
        // This is simplified - in practice we'd identify SNVs through alignment to reference

        Ok(())
    }

    /// Calculate weighted similarity between this signature and another
    pub fn similarity(
        &self,
        other: &MultiResolutionSignature,
        level: Option<ResolutionLevel>,
    ) -> f64 {
        match level {
            Some(ResolutionLevel::Macro) => self
                .macro_signature
                .jaccard_similarity(&other.macro_signature),
            Some(ResolutionLevel::Meso) => self
                .meso_signature
                .jaccard_similarity(&other.meso_signature),
            Some(ResolutionLevel::Micro) => self.micro_signature.similarity(&other.micro_signature),
            None => {
                // Weighted combination of all levels
                let macro_sim = self
                    .macro_signature
                    .jaccard_similarity(&other.macro_signature);
                let meso_sim = self
                    .meso_signature
                    .jaccard_similarity(&other.meso_signature);
                let micro_sim = self.micro_signature.similarity(&other.micro_signature);

                // Get weights or use defaults
                let macro_weight = self.weights.get(&ResolutionLevel::Macro).unwrap_or(&0.3);
                let meso_weight = self.weights.get(&ResolutionLevel::Meso).unwrap_or(&0.4);
                let micro_weight = self.weights.get(&ResolutionLevel::Micro).unwrap_or(&0.3);

                // Calculate weighted average
                let total_weight = macro_weight + meso_weight + micro_weight;
                if total_weight > 0.0 {
                    (macro_sim * macro_weight + meso_sim * meso_weight + micro_sim * micro_weight)
                        / total_weight
                } else {
                    (macro_sim + meso_sim + micro_sim) / 3.0
                }
            }
        }
    }

    /// Calculate information content (discriminatory power) of this signature
    pub fn calculate_weights(&mut self, references: &[&MultiResolutionSignature]) {
        // Implementation of information-theoretic weighting
        // This is a simplified approach - full implementation would be more sophisticated

        let mut macro_differences = Vec::new();
        let mut meso_differences = Vec::new();
        let mut micro_differences = Vec::new();

        for reference in references {
            if reference.taxon_id == self.taxon_id {
                continue; // Skip self-comparison
            }

            // Calculate differences at each level
            let macro_diff = 1.0
                - self
                    .macro_signature
                    .jaccard_similarity(&reference.macro_signature);
            let meso_diff = 1.0
                - self
                    .meso_signature
                    .jaccard_similarity(&reference.meso_signature);
            let micro_diff = 1.0 - self.micro_signature.similarity(&reference.micro_signature);

            macro_differences.push(macro_diff);
            meso_differences.push(meso_diff);
            micro_differences.push(micro_diff);
        }

        // Calculate average differences (higher = more discriminatory)
        let macro_avg = if !macro_differences.is_empty() {
            macro_differences.iter().sum::<f64>() / macro_differences.len() as f64
        } else {
            0.3 // Default
        };

        let meso_avg = if !meso_differences.is_empty() {
            meso_differences.iter().sum::<f64>() / meso_differences.len() as f64
        } else {
            0.4 // Default
        };

        let micro_avg = if !micro_differences.is_empty() {
            micro_differences.iter().sum::<f64>() / micro_differences.len() as f64
        } else {
            0.3 // Default
        };

        // Normalize weights to sum to 1.0
        let total = macro_avg + meso_avg + micro_avg;
        if total > 0.0 {
            self.weights
                .insert(ResolutionLevel::Macro, macro_avg / total);
            self.weights.insert(ResolutionLevel::Meso, meso_avg / total);
            self.weights
                .insert(ResolutionLevel::Micro, micro_avg / total);
        } else {
            // Default weights if no discrimination power
            self.weights.insert(ResolutionLevel::Macro, 0.3);
            self.weights.insert(ResolutionLevel::Meso, 0.4);
            self.weights.insert(ResolutionLevel::Micro, 0.3);
        }
    }
}

/// Builder for creating multi-resolution signatures from genomes
pub struct SignatureBuilder {
    /// K-mer size for macro-level (species)
    macro_k: usize,

    /// K-mer size for meso-level (strain group)
    meso_k: usize,

    /// Sketch size for MinHash signatures
    sketch_size: usize,

    /// Number of threads for parallel processing
    threads: usize,
}

impl SignatureBuilder {
    /// Create a new signature builder
    pub fn new(
        macro_k: usize,
        meso_k: usize,
        sketch_size: usize,
        threads: usize,
    ) -> Result<Self, SignatureError> {
        if macro_k < 3 || macro_k > 31 || meso_k < 3 || meso_k > 31 {
            return Err(SignatureError::InvalidKmerSize(
                if macro_k < 3 || macro_k > 31 {
                    macro_k
                } else {
                    meso_k
                },
            ));
        }

        Ok(SignatureBuilder {
            macro_k,
            meso_k,
            sketch_size,
            threads: threads.max(1),
        })
    }

    /// Build a signature from a FASTA file
    pub fn build_from_file<P: AsRef<Path>>(
        &self,
        file_path: P,
        taxon_id: &str,
        lineage: Vec<String>,
    ) -> Result<MultiResolutionSignature, SignatureError> {
        // Create empty signature
        let mut signature = MultiResolutionSignature::new(
            taxon_id,
            lineage,
            self.macro_k,
            self.meso_k,
            self.sketch_size,
            self.sketch_size,
        )?;

        // Parse FASTA file
        let mut reader = parse_fastx_file(file_path).map_err(|e| {
            SignatureError::IoError(std::io::Error::new(std::io::ErrorKind::Other, e))
        })?;

        // Process each sequence in the file
        while let Some(record) = reader.next() {
            let record = record.map_err(|e| {
                SignatureError::IoError(std::io::Error::new(std::io::ErrorKind::Other, e))
            })?;

            // Normalize sequence to uppercase and filter out non-ACGT
            let seq = record.seq();
            signature.add_sequence(&seq)?;
        }

        Ok(signature)
    }

    /// Build multiple signatures in parallel
    pub fn build_batch<P: AsRef<Path> + Sync>(
        &self,
        files: Vec<(P, String, Vec<String>)>, // (path, taxon_id, lineage)
    ) -> Result<Vec<MultiResolutionSignature>, SignatureError> {
        // Process files in parallel using rayon
        let results: Vec<Result<MultiResolutionSignature, SignatureError>> = files
            .par_iter()
            .map(|(path, taxon_id, lineage)| self.build_from_file(path, taxon_id, lineage.clone()))
            .collect();

        // Collect results, propagating any errors
        let mut signatures = Vec::with_capacity(results.len());
        for result in results {
            signatures.push(result?);
        }

        // Calculate weights based on discriminatory power
        for i in 0..signatures.len() {
            let (signature, others) = get_signature_and_others(&mut signatures, i);
            let others_refs: Vec<&MultiResolutionSignature> = others.iter().cloned().collect();
            signature.calculate_weights(&others_refs);
        }

        Ok(signatures)
    }
}

/// Helper function to get a mutable reference to a signature and slice of others
fn get_signature_and_others<'a>(
    signatures: &'a mut Vec<MultiResolutionSignature>,
    index: usize,
) -> (
    &'a mut MultiResolutionSignature,
    Vec<&'a MultiResolutionSignature>,
) {
    let mut others = Vec::with_capacity(signatures.len() - 1);

    // Split the vector to get a mutable reference to the target and references to others
    let (left, right) = signatures.split_at_mut(index);
    let (signature, right_rest) = right.split_first_mut().unwrap();

    // Collect references to other signatures
    for sig in left.iter() {
        others.push(sig);
    }
    for sig in right_rest.iter() {
        others.push(sig);
    }

    (signature, others)
}
