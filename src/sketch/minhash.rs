//! MinHash sketching implementation.
//!
//! MinHash is a technique to estimate the Jaccard similarity between sets
//! (in this case, sets of k-mers derived from sequences) quickly using
//! fixed-size sketches (signatures).

use crate::sketch::signature::Signature;
use crate::sketch::Sketcher; // Implement the common Sketcher trait

use anyhow::{anyhow, Result};
use bio::io::fasta::Record as SequenceRecord2;
use needletail::parser::SequenceRecord;
use needletail::FastxReader;
use needletail::Sequence;
use std::collections::hash_map::DefaultHasher; // Simple default hasher
use std::hash::{Hash, Hasher};
// Consider using more robust hashing like xxHash or MurmurHash3 via crates
// e.g., use fasthash::xx;

/// Structure for creating MinHash sketches.
#[derive(Debug, Clone)]
pub struct MinHashSketcher {
    num_hashes: usize, // Number of hash values in the sketch (sketch size)
    kmer_size: usize,  // K-mer size to use
                       // TODO: Potentially add seeds or precomputed hash functions if not using a single hasher.
                       // seeds: Vec<u64>,
}

impl MinHashSketcher {
    /// Creates a new MinHashSketcher.
    ///
    /// # Arguments
    /// * `num_hashes` - The desired sketch size (number of hash values).
    /// * `kmer_size` - The k-mer length.
    pub fn new(num_hashes: usize, kmer_size: usize) -> Result<Self> {
        if num_hashes == 0 {
            return Err(anyhow!(
                "Number of hashes (sketch size) must be greater than 0."
            ));
        }
        if kmer_size == 0 {
            return Err(anyhow!("K-mer size must be greater than 0."));
        }
        Ok(MinHashSketcher {
            num_hashes,
            kmer_size,
        })
    }

    /// Calculates a hash value for a k-mer.
    /// Uses the standard library's DefaultHasher for simplicity.
    // TODO: Replace with a more robust/faster hash function if needed.
    fn hash_kmer(&self, kmer: &[u8]) -> u64 {
        let mut hasher = DefaultHasher::new();
        kmer.hash(&mut hasher);
        hasher.finish()
    }

    // TODO: If using multiple hash functions (more standard MinHash):
    // fn hash_kmer_seeded(&self, kmer: &[u8], seed: u64) -> u64 {
    //     let mut hasher = // Initialize hasher with seed, e.g., xx::Hasher64::with_seed(seed);
    //     kmer.hash(&mut hasher);
    //     hasher.finish()
    // }
}

impl Sketcher for MinHashSketcher {
    /// Creates a MinHash signature for a single sequence record.
    /// This implementation uses the "bottom-k" approach with a single hash function
    /// for simplicity, which is technically equivalent to MinHash under certain assumptions.
    /// A more standard implementation would use `num_hashes` different hash functions.
    fn sketch_sequence(&self, record: &SequenceRecord) -> Result<Signature> {
        let mut signature = Signature::new("minhash".to_string(), self.kmer_size, self.num_hashes);
        signature.name = Some(String::from_utf8_lossy(record.id()).into_owned());
        // signature.filename = ... // Can be set later if sketching from a file context

        let seq = record.sequence();
        let mut kmer_hashes = Vec::new();

        // 1. Generate all canonical k-mers and hash them
        // TODO: Use the actual CanonicalKmerIter from bio::kmers when implemented correctly.
        // Manual canonicalization for now:
        for i in 0..=(seq.len().saturating_sub(self.kmer_size)) {
            let kmer = &seq[i..i + self.kmer_size];
            if kmer.iter().any(|&b| !crate::bio::is_valid_base(b)) {
                continue; // Skip k-mers with invalid bases
            }
            let rc = crate::bio::reverse_complement(kmer);
            let canonical_kmer = if kmer < &rc[..] { kmer } else { &rc[..] };
            kmer_hashes.push(self.hash_kmer(canonical_kmer));
        }

        // 2. Keep the smallest `num_hashes` unique hash values (Bottom-k sketch)
        kmer_hashes.sort_unstable(); // Sort hashes
        kmer_hashes.dedup(); // Remove duplicates

        // Take the smallest `num_hashes` values
        signature.hashes = kmer_hashes.into_iter().take(self.num_hashes).collect();

        // Ensure the signature hash count matches num_hashes if fewer unique hashes were found
        // This is important for consistent Jaccard estimation later.
        // Standard MinHash often pads with MAX_U64, but just storing fewer is also common.
        // Let's store fewer for now, but update num_hashes in the signature itself might be better.
        // signature.num_hashes = signature.hashes.len(); // Alternative: adjust signature metadata

        Ok(signature)
    }

    // Override sketch_sequences for potential optimization if needed,
    // otherwise the default implementation from the trait is used.
}
