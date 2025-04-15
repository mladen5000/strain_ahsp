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
use signature::{KmerSignature, KmerSignatureBuilder};

// Re-export key structures or functions if needed
// pub use minhash::MinHashSketcher;
// pub use adaptive::AdaptiveSketcher;

use crate::sketch::signature::Signature; // Use our own Signature structure
use anyhow::{anyhow, Result};
use needletail::parser::SequenceRecord;
use needletail::{parse_fastx_file, Sequence};
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
    pub fn build_from_file<P: AsRef<Path>>(
        &self,
        file_path: P,
        taxon_id: &str,
        lineage: Vec<String>,
    ) -> Result<MultiResolutionSignature> {
        let mut reader = parse_fastx_file(file_path.as_ref())?;
        let mut multi_sig = MultiResolutionSignature::new(taxon_id.to_string(), lineage);

        // Calculate k-mer sizes for each level
        let k_step = (self.kmer_size - self.min_kmer_size) as f32 / (self.levels - 1) as f32;

        // Create signatures for each resolution level
        for level in 0..self.levels {
            let level_k = (self.kmer_size as f32 - (level as f32 * k_step)).round() as usize;
            let level_sketch_size = self.sketch_size / (1 << level); // Decrease sketch size for finer resolutions

            let mut level_sig =
                KmerSignatureBuilder::new(level_k, "DNA", "minhash", level_sketch_size, 0);
            let name = format!("level_{}", level);
            let built_sig = level_sig.source(file_path.as_ref()).name(&name).build();
            multi_sig.add_level(built_sig);
        }

        Ok(multi_sig)
    }

    /// Builds multiple signatures from a batch of FASTA/FASTQ files.
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
