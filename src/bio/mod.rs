//! Bioinformatics utilities module.
//!
//! This module likely groups together sub-modules related to
//! biological sequence processing, k-mer manipulation, and
//! potentially signature/sketch generation specific details.

// Declare sub-modules within the 'bio' directory
pub mod kmers;
pub mod signature; // Module for handling sequence signatures (e.g., from sketching)
pub mod taxonomy;

pub use kmers::KmerExtractor;
pub use taxonomy::{TaxonomicLevel, TaxonomicLineage};

// Re-export important items from sub-modules if desired
// pub use kmers::Kmer;
// pub use signature::Signature;

/// General bioinformatics constants or utility functions can be placed here.
pub const CANONICAL_BASES: &[u8] = b"ACGT";

/// Checks if a byte represents a valid DNA base (A, C, G, T).
/// Case-insensitive.
pub fn is_valid_base(base: u8) -> bool {
    matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
}

/// Calculates the reverse complement of a DNA sequence.
/// Handles IUPAC codes partially (N -> N). Others might become N or cause errors.
pub fn reverse_complement(dna: &[u8]) -> Vec<u8> {
    dna.iter()
        .rev()
        .map(|&base| match base.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'N' => b'N', // Keep N as N
            _ => b'N',    // Or return an error / handle other IUPAC codes if needed
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_valid_base_standard() {
        assert!(is_valid_base(b'A'));
        assert!(is_valid_base(b'c'));
        assert!(is_valid_base(b'G'));
        assert!(is_valid_base(b't'));
    }

    #[test]
    fn test_is_valid_base_invalid() {
        assert!(!is_valid_base(b'N')); // N is often handled separately
        assert!(!is_valid_base(b'X'));
        assert!(!is_valid_base(b' '));
    }

    #[test]
    fn test_reverse_complement_simple() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b"GATTACA"), b"TGTAATC");
    }

    #[test]
    fn test_reverse_complement_mixed_case() {
        assert_eq!(reverse_complement(b"aCgT"), b"ACGT");
    }

    #[test]
    fn test_reverse_complement_with_n() {
        assert_eq!(reverse_complement(b"ANT"), b"ANT");
    }

    #[test]
    fn test_reverse_complement_empty() {
        assert_eq!(reverse_complement(b""), b"");
    }
}
