//! Module for sequence sketching algorithms.
//!
//! Sketching techniques like MinHash or KHF (K-mer Hashing Frameworks)
//! allow for fast comparison and clustering of large sequences or datasets
//! by creating compressed representations (signatures or sketches).

pub mod adaptive;
pub mod minhash; // MinHash implementation // Potentially adaptive MinHash or other adaptive sketching
pub mod signature;

pub use adaptive::AdaptiveClassifier;
pub use signature::MultiResolutionSignature;

// Re-export key structures or functions if needed
// pub use minhash::MinHashSketcher;
// pub use adaptive::AdaptiveSketcher;

use crate::sketch::signature::Signature; // Use our own Signature structure
use anyhow::{anyhow, Result};
use needletail::parse_fastx_file;
use needletail::parser::SequenceRecord;
use std::path::Path;

/// Trait defining common operations for sequence sketchers.
pub trait Sketcher {
    /// Creates a signature (sketch) for a single sequence record.
    ///
    /// # Arguments
    /// * `record` - The sequence record to sketch.
    ///
    /// # Returns
    /// * `Result<Signature>` - The generated signature or an error.
    fn sketch_sequence(&self, record: &SequenceRecord) -> Result<Signature>;

    /// Creates signatures for multiple sequence records.
    /// Could be implemented more efficiently than calling `sketch_sequence` repeatedly.
    ///
    /// # Arguments
    /// * `records` - An iterator over sequence records.
    ///
    /// # Returns
    /// * `Result<Vec<Signature>>` - A vector of generated signatures or an error.
    fn sketch_sequences<'a>(
        &self,
        records: impl Iterator<Item = Result<SequenceRecord<'a>>>,
    ) -> Result<Vec<Signature>> {
        // Default implementation: iterate and call sketch_sequence
        let mut signatures = Vec::new();
        for record_result in records {
            let record = record_result?;
            signatures.push(self.sketch_sequence(&record)?);
        }
        Ok(signatures)
    }

    // TODO: Add methods for sketching entire files or combining sketches if applicable.
}

/// Builder for creating genomic signatures from sequence data.
pub struct SignatureBuilder {
    pub kmer_size: u8,
    pub sketch_size: usize,
    pub min_kmer_size: u8,
    pub levels: u8,
}

impl SignatureBuilder {
    /// Creates a new SignatureBuilder with the specified parameters.
    ///
    /// # Arguments
    ///
    /// * `kmer_size` - The size of k-mers to use (e.g., 31)
    /// * `min_kmer_size` - The minimum k-mer size for adaptive resolution (e.g., 21)
    /// * `sketch_size` - The size of each sketch (e.g., 1000)
    /// * `levels` - The number of hierarchical levels to generate (e.g., 8)
    ///
    /// # Returns
    ///
    /// A Result containing the new SignatureBuilder or an error
    pub fn new(kmer_size: u8, min_kmer_size: u8, sketch_size: usize, levels: u8) -> Result<Self> {
        if kmer_size < min_kmer_size {
            return Err(anyhow!("kmer_size must be >= min_kmer_size"));
        }
        if kmer_size > 63 {
            return Err(anyhow!("kmer_size must be <= 63"));
        }
        if levels == 0 {
            return Err(anyhow!("levels must be > 0"));
        }

        Ok(Self {
            kmer_size,
            sketch_size,
            min_kmer_size,
            levels,
        })
    }

    /// Builds a signature from a FASTA/FASTQ file.
    ///
    /// # Arguments
    ///
    /// * `file_path` - Path to the input sequence file
    /// * `taxon_id` - Identifier for the sequence
    /// * `lineage` - Taxonomic lineage (if available)
    ///
    /// # Returns
    ///
    /// A Result containing the new MultiResolutionSignature or an error
    pub fn build_from_file<P: AsRef<Path>>(
        &self,
        file_path: P,
        taxon_id: &str,
        lineage: Vec<String>,
    ) -> Result<MultiResolutionSignature> {
        // Implementation would typically:
        // 1. Parse the file using needletail
        // 2. Extract k-mers and compute their hashes
        // 3. Build the hierarchical sketches
        // 4. Return a MultiResolutionSignature

        // Placeholder implementation
        let mut reader = parse_fastx_file(file_path)?;

        // This is just a placeholder - actual implementation would process the sequence
        Ok(MultiResolutionSignature::new(taxon_id.to_string(), lineage))
    }

    /// Builds multiple signatures from a batch of FASTA/FASTQ files.
    ///
    /// # Arguments
    ///
    /// * `files` - Vector of tuples containing (file_path, taxon_id, lineage)
    ///
    /// # Returns
    ///
    /// A Result containing a vector of MultiResolutionSignature or an error
    pub fn build_batch<P: AsRef<Path>>(
        &self,
        files: Vec<(P, String, Vec<String>)>,
    ) -> Result<Vec<MultiResolutionSignature>> {
        let mut signatures = Vec::with_capacity(files.len());

        for (file_path, taxon_id, lineage) in files {
            let signature = self.build_from_file(file_path, &taxon_id, lineage)?;
            signatures.push(signature);
        }

        Ok(signatures)
    }
}

// Re-exported above

#[cfg(test)]
mod tests {
    // Add tests for any functions or constants defined directly in this mod.rs file.
    // Tests for specific sketchers should go into their respective modules (minhash.rs, adaptive.rs).
}
