//! MinHash sketching implementation.
//!
//! MinHash is a technique to estimate the Jaccard similarity between sets
//! (in this case, sets of k-mers derived from sequences) quickly using
//! fixed-size sketches (signatures).

use crate::sketch::signature::Signature;
use crate::sketch::Sketcher; // Implement the common Sketcher trait

use anyhow::{anyhow, Result};
use needletail::parser::SequenceRecord;
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

#[cfg(test)]
mod tests {
    use super::*;
    use needletail::parser::FastxReader;
    use needletail::parser::SequenceRecord;
    use std::io::Cursor;

    // Helper to create a SequenceRecord from a string
    fn seq_rec(id: &'static str, seq: &'static str) -> SequenceRecord<'static> {
        let data = format!(">{}\n{}\n", id, seq);
        let cursor = Cursor::new(data.into_bytes());
        let mut reader = FastxReader::new(cursor);
        reader
            .next()
            .unwrap()
            .expect("Failed to parse test sequence")
    }

    #[test]
    fn test_minhash_sketcher_new() {
        let sketcher = MinHashSketcher::new(100, 21).unwrap();
        assert_eq!(sketcher.num_hashes, 100);
        assert_eq!(sketcher.kmer_size, 21);
    }

    #[test]
    fn test_minhash_sketcher_new_invalid() {
        assert!(MinHashSketcher::new(0, 21).is_err());
        assert!(MinHashSketcher::new(100, 0).is_err());
    }

    #[test]
    fn test_sketch_sequence_basic() {
        // Use a small k and num_hashes for predictability (though hashing makes it tricky)
        let sketcher = MinHashSketcher::new(5, 3).unwrap(); // num_hashes=5, k=3
        let record = seq_rec("test_seq", "ACGTTACGT");
        // 3-mers: ACG, CGT, GTT, TTA, TAC, ACG, CGT
        // Canonical: ACG, ACG, AAC(rev GTT), TAA(rev TTA), TAC, ACG, ACG
        // Hashes (example, depends on hasher): h(ACG), h(AAC), h(TAA), h(TAC)
        // Assume hashes are unique and sorted: [h(AAC), h(ACG), h(TAA), h(TAC)]
        // Sketch (num_hashes=5, but only 4 unique): [h(AAC), h(ACG), h(TAA), h(TAC)]

        let signature_res = sketcher.sketch_sequence(&record);
        assert!(signature_res.is_ok());
        let signature = signature_res.unwrap();

        assert_eq!(signature.algorithm, "minhash");
        assert_eq!(signature.kmer_size, 3);
        assert_eq!(signature.num_hashes, 5); // Original requested size
        assert_eq!(signature.name.unwrap(), "test_seq");

        // Verify the resulting hashes (bottom-k of unique canonical k-mer hashes)
        // This requires knowing the hash function precisely. We'll check the count.
        let mut expected_hashes = vec![
            sketcher.hash_kmer(b"AAC"), // rev(GTT)
            sketcher.hash_kmer(b"ACG"), // ACG and rev(CGT)
            sketcher.hash_kmer(b"TAA"), // rev(TTA)
            sketcher.hash_kmer(b"TAC"), // TAC and rev(GTA) - assuming TAC > GTA
        ];
        expected_hashes.sort_unstable();
        expected_hashes.dedup(); // Should already be unique

        let mut actual_hashes = signature.hashes;
        actual_hashes.sort_unstable(); // Sort actual hashes for comparison

        assert_eq!(actual_hashes.len(), 4); // Should contain only the unique hashes found
        assert_eq!(actual_hashes, expected_hashes);
    }

    #[test]
    fn test_sketch_sequence_k_too_large() {
        let sketcher = MinHashSketcher::new(10, 5).unwrap();
        let record = seq_rec("short", "ACGT"); // k=5, seq len=4
        let signature = sketcher.sketch_sequence(&record).unwrap();
        assert!(signature.hashes.is_empty());
    }

    #[test]
    fn test_sketch_sequence_with_n() {
        let sketcher = MinHashSketcher::new(10, 3).unwrap();
        let record = seq_rec("with_n", "ACNGT");
        // 3-mers: ACN, CNG, NGT -> all skipped
        let signature = sketcher.sketch_sequence(&record).unwrap();
        assert!(signature.hashes.is_empty());
    }

    #[test]
    fn test_sketch_sequence_fewer_hashes_than_num() {
        // Only 1 unique canonical 3-mer (AAA)
        let sketcher = MinHashSketcher::new(10, 3).unwrap();
        let record = seq_rec("repeats", "AAAAAAA");
        let signature = sketcher.sketch_sequence(&record).unwrap();

        let expected_hash = sketcher.hash_kmer(b"AAA");
        assert_eq!(signature.hashes.len(), 1);
        assert_eq!(signature.hashes[0], expected_hash);
        assert_eq!(signature.num_hashes, 10); // Metadata should still reflect requested size
    }
}
