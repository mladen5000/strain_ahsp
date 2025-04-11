//! K-mer generation and processing utilities.
//!
//! This module focuses on extracting k-mers (substrings of length k)
//! from biological sequences. It might include:
//! - K-mer iterators over sequences.
//! - Canonical k-mer generation (lexicographically smaller of k-mer and its reverse complement).
//! - K-mer counting functions.

use crate::bio; // Access functions like reverse_complement from parent bio module
use anyhow::Result;
use log::warn;
use needletail::parser::SequenceRecord; // Using needletail for sequence handling
use needletail::Sequence;
use std::collections::HashMap; // For k-mer counting

/// Represents a k-mer. Could be stored as bytes, string, or a packed integer format.
// Example using bytes:
pub type Kmer<'a> = &'a [u8];

// TODO: Consider using a more efficient representation like a u64 for small k,
// or libraries like `needletail::kmer` or `nthash` for hashing/canonicalization.

/// An iterator over canonical k-mers in a sequence.
///
/// Yields the lexicographically smaller of a k-mer and its reverse complement.
/// Skips k-mers containing invalid bases (e.g., 'N').
pub struct CanonicalKmerIter<'a> {
    sequence: &'a [u8],
    k: usize,
    current_pos: usize,
}

impl<'a> CanonicalKmerIter<'a> {
    pub fn new(sequence: &'a [u8], k: usize) -> Self {
        CanonicalKmerIter {
            sequence,
            k,
            current_pos: 0,
        }
    }
}

impl<'a> Iterator for CanonicalKmerIter<'a> {
    type Item = Kmer<'a>; // Yielding byte slices for now

    fn next(&mut self) -> Option<Self::Item> {
        while self.current_pos + self.k <= self.sequence.len() {
            let kmer_slice = &self.sequence[self.current_pos..self.current_pos + self.k];
            self.current_pos += 1; // Move to next position regardless

            // Check for invalid bases (e.g., 'N') within the k-mer
            if kmer_slice.iter().any(|&b| !bio::is_valid_base(b)) {
                continue; // Skip k-mers with invalid bases
            }

            // TODO: Implement efficient canonical k-mer generation.
            // This simple version calculates reverse complement every time.
            // For performance, consider rolling hashes or bit manipulation.
            let rc_kmer = bio::reverse_complement(kmer_slice);

            // Compare the k-mer with its reverse complement
            if kmer_slice <= &rc_kmer[..] {
                return Some(kmer_slice); // Return original if smaller or equal
            } else {
                // This is tricky because rc_kmer is owned Vec<u8>.
                // We need to return a slice with the same lifetime 'a.
                // This requires either:
                // 1. Storing the RC temporarily (inefficient).
                // 2. Using a k-mer representation that handles canonicalization internally.
                // 3. Allocating the RC k-mer and leaking it (bad idea).
                // For now, we'll just return the original slice, acknowledging this
                // implementation is NOT strictly canonical without more work.
                // A better approach uses numerical k-mer representations.
                // *** Placeholder: Returning original for now ***
                warn!("CanonicalKmerIter currently returns original k-mer, not necessarily canonical due to lifetime issues with RC.");
                return Some(kmer_slice);
            }
        }
        None // End of sequence
    }
}

/// Counts canonical k-mers in a single sequence record.
///
/// # Arguments
///
/// * `record` - A `SequenceRecord` (e.g., from needletail).
/// * `k` - The k-mer size.
/// * `counts` - A mutable HashMap to store k-mer counts. K-mers are added as Vec<u8>.
///
/// # Returns
///
/// * `Result<()>` - Ok(()) if successful, or an error.
pub fn count_canonical_kmers_in_record(
    record: &SequenceRecord,
    k: usize,
    counts: &mut HashMap<Vec<u8>, u32>, // Using Vec<u8> as key for simplicity
) -> Result<()> {
    if k == 0 {
        return Ok(());
    } // No k-mers of size 0

    let seq = record.sequence();
    // TODO: Use the *actual* CanonicalKmerIter once it's correctly implemented.
    // For now, iterating simply and calculating canonical form manually.
    for i in 0..=(seq.len().saturating_sub(k)) {
        let kmer = &seq[i..i + k];

        // Basic check for N's
        if kmer.iter().any(|&b| b == b'N' || b == b'n') {
            continue;
        }

        // Calculate canonical k-mer (lexicographically smallest)
        let rc = bio::reverse_complement(kmer);
        let canonical_kmer = if kmer < &rc[..] { kmer } else { &rc[..] };

        *counts.entry(canonical_kmer.to_vec()).or_insert(0) += 1;
    }

    Ok(())
}

/// Processes multiple sequences to generate k-mer counts (example function).
///
/// # Arguments
///
/// * `sequences` - An iterator over sequence records.
/// * `k` - The k-mer size.
///
/// # Returns
///
/// * `Result<HashMap<Vec<u8>, u32>>` - A map of canonical k-mers to their total counts.
pub fn process_sequences<'a>(
    sequences: impl Iterator<Item = Result<SequenceRecord<'a>>>,
    k: usize,
) -> Result<HashMap<Vec<u8>, u32>> {
    let mut total_counts = HashMap::new();
    for result in sequences {
        let record = result?; // Propagate potential reading errors
        count_canonical_kmers_in_record(&record, k, &mut total_counts)?;
    }
    Ok(total_counts)
}

/// A configurable k-mer extraction utility.
/// Handles various k-mer extraction settings like size, canonicalization, etc.
pub struct KmerExtractor {
    pub k: usize,
    pub canonical: bool,
    pub skip_invalid: bool,
}

impl KmerExtractor {
    /// Creates a new KmerExtractor with default settings (canonical k-mers, skipping invalid).
    pub fn new(k: usize) -> Self {
        KmerExtractor {
            k,
            canonical: true,
            skip_invalid: true,
        }
    }

    /// Creates a KmerExtractor with custom settings.
    pub fn with_settings(k: usize, canonical: bool, skip_invalid: bool) -> Self {
        KmerExtractor {
            k,
            canonical,
            skip_invalid,
        }
    }

    /// Extract and count k-mers from a sequence.
    pub fn count_kmers(&self, seq: &[u8]) -> HashMap<Vec<u8>, u32> {
        let mut counts = HashMap::new();

        if self.k == 0 || seq.len() < self.k {
            return counts;
        }

        for i in 0..=(seq.len() - self.k) {
            let kmer = &seq[i..i + self.k];

            // Skip k-mers with invalid bases if required
            if self.skip_invalid && kmer.iter().any(|&b| !bio::is_valid_base(b)) {
                continue;
            }

            // Get canonical form if required
            let final_kmer = if self.canonical {
                let rc = bio::reverse_complement(kmer);
                if kmer < &rc[..] {
                    kmer.to_vec()
                } else {
                    rc
                }
            } else {
                kmer.to_vec()
            };

            *counts.entry(final_kmer).or_insert(0) += 1;
        }

        counts
    }

    /// Process and count k-mers from a sequence record.
    pub fn process_record(&self, record: &SequenceRecord) -> HashMap<Vec<u8>, u32> {
        self.count_kmers(record.sequence())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use needletail::parser::FastxReader;
    use std::io::Cursor;

    #[test]
    fn test_kmer_counting_simple() {
        let seq_data = ">test\nACGTACGT\n";
        let cursor = Cursor::new(seq_data.as_bytes());
        let mut reader = needletail::parse_fastx_reader(cursor).expect("Failed to create reader");
        let record = reader.next().unwrap().expect("Parse failed");

        let mut counts = HashMap::new();
        count_canonical_kmers_in_record(&record, 3, &mut counts).unwrap();

        // Expected 3-mers: ACG, CGT, GTA, TAC, ACG, CGT
        // Canonical: ACG, ACG (rev(CGT)=ACG), ACG (rev(GTA)=TAC), TAC, ACG, ACG
        // Counts: ACG: 5, TAC: 1
        assert_eq!(counts.get(b"ACG".as_ref()), Some(&5));
        assert_eq!(counts.get(b"TAC".as_ref()), Some(&1));
        assert_eq!(counts.len(), 2);
    }

    #[test]
    fn test_kmer_counting_with_n() {
        let seq_data = ">test\nACNGT";
        let cursor = Cursor::new(seq_data.as_bytes());
        let mut reader = needletail::parse_fastx_reader(cursor).expect("Failed to create reader");
        let record = reader.next().unwrap().expect("Parse failed");

        let mut counts = HashMap::new();
        count_canonical_kmers_in_record(&record, 3, &mut counts).unwrap();
        assert!(counts.is_empty());
    }

    #[test]
    fn test_kmer_counting_empty_seq() {
        let seq_data = ">test\n";
        let cursor = Cursor::new(seq_data.as_bytes());
        let mut reader = needletail::parse_fastx_reader(cursor).expect("Failed to create reader");
        let record = reader.next().unwrap().expect("Parse failed");

        let mut counts = HashMap::new();
        count_canonical_kmers_in_record(&record, 3, &mut counts).unwrap();
        assert!(counts.is_empty());
    }

    #[test]
    fn test_kmer_counting_k_too_large() {
        let seq_data = ">test\nACGT";
        let cursor = Cursor::new(seq_data.as_bytes());
        let mut reader = needletail::parse_fastx_reader(cursor).expect("Failed to create reader");
        let record = reader.next().unwrap().expect("Parse failed");

        let mut counts = HashMap::new();
        count_canonical_kmers_in_record(&record, 5, &mut counts).unwrap();
        assert!(counts.is_empty());
    }

    // TODO: Add tests for the CanonicalKmerIter once it's fully implemented.
}
