use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use thiserror::Error;

use crate::sketch::signature::{MultiResolutionSignature, ResolutionLevel};

#[derive(Error, Debug)]
pub enum ClassificationError {
    #[error("No reference signatures available")]
    NoReferences,
    #[error("Insufficient coverage for classification")]
    InsufficientCoverage,
}

/// Taxonomic levels from domain to strain
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum TaxonomicLevel {
    Domain,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species,
    StrainGroup,
    Strain,
    Unknown,
}

impl TaxonomicLevel {
    /// Get parent taxonomic level
    pub fn parent(&self) -> Self {
        match self {
            TaxonomicLevel::Domain => TaxonomicLevel::Unknown,
            TaxonomicLevel::Phylum => TaxonomicLevel::Domain,
            TaxonomicLevel::Class => TaxonomicLevel::Phylum,
            TaxonomicLevel::Order => TaxonomicLevel::Class,
            TaxonomicLevel::Family => TaxonomicLevel::Order,
            TaxonomicLevel::Genus => TaxonomicLevel::Family,
            TaxonomicLevel::Species => TaxonomicLevel::Genus,
            TaxonomicLevel::StrainGroup => TaxonomicLevel::Species,
            TaxonomicLevel::Strain => TaxonomicLevel::StrainGroup,
            TaxonomicLevel::Unknown => TaxonomicLevel::Unknown,
        }
    }

    /// Get index in lineage array (0-based)
    pub fn lineage_index(&self) -> Option<usize> {
        match self {
            TaxonomicLevel::Domain => Some(0),
            TaxonomicLevel::Phylum => Some(1),
            TaxonomicLevel::Class => Some(2),
            TaxonomicLevel::Order => Some(3),
            TaxonomicLevel::Family => Some(4),
            TaxonomicLevel::Genus => Some(5),
            TaxonomicLevel::Species => Some(6),
            TaxonomicLevel::StrainGroup => Some(7),
            TaxonomicLevel::Strain => Some(8),
            TaxonomicLevel::Unknown => None,
        }
    }
}

/// Confidence thresholds for different taxonomic levels
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConfidenceThresholds {
    pub thresholds: HashMap<TaxonomicLevel, f64>,
}

impl Default for ConfidenceThresholds {
    fn default() -> Self {
        let mut thresholds = HashMap::new();
        thresholds.insert(TaxonomicLevel::Domain, 0.95);
        thresholds.insert(TaxonomicLevel::Phylum, 0.92);
        thresholds.insert(TaxonomicLevel::Class, 0.90);
        thresholds.insert(TaxonomicLevel::Order, 0.87);
        thresholds.insert(TaxonomicLevel::Family, 0.85);
        thresholds.insert(TaxonomicLevel::Genus, 0.80);
        thresholds.insert(TaxonomicLevel::Species, 0.75);
        thresholds.insert(TaxonomicLevel::StrainGroup, 0.70);
        thresholds.insert(TaxonomicLevel::Strain, 0.65);

        ConfidenceThresholds { thresholds }
    }
}

/// Classification result
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Classification {
    /// Taxonomic ID at the assigned level
    pub taxon_id: String,

    /// Full lineage from domain to strain (when available)
    pub lineage: Vec<String>,

    /// Taxonomic level of the classification
    pub level: TaxonomicLevel,

    /// Classification confidence
    pub confidence: f64,

    /// Best matching reference signature
    pub best_match: String,

    /// Similarity scores at each resolution level
    pub similarity_scores: HashMap<ResolutionLevel, f64>,
}

/// Adaptive resolution classifier
pub struct AdaptiveClassifier {
    /// Reference signatures
    pub references: Vec<MultiResolutionSignature>,

    /// Confidence thresholds
    pub thresholds: ConfidenceThresholds,

    /// Mapping of reference IDs to indices
    reference_index: HashMap<String, usize>,

    /// Minimum coverage required for classification
    min_coverage: usize,
}

impl AdaptiveClassifier {
    /// Create a new adaptive classifier
    pub fn new(
        references: Vec<MultiResolutionSignature>,
        thresholds: Option<ConfidenceThresholds>,
        min_coverage: Option<usize>,
    ) -> Result<Self, ClassificationError> {
        if references.is_empty() {
            return Err(ClassificationError::NoReferences);
        }

        // Build reference index
        let mut reference_index = HashMap::new();
        for (i, ref_sig) in references.iter().enumerate() {
            reference_index.insert(ref_sig.taxon_id.clone(), i);
        }

        Ok(AdaptiveClassifier {
            references,
            thresholds: thresholds.unwrap_or_default(),
            reference_index,
            min_coverage: min_coverage.unwrap_or(100),
        })
    }

    /// Classify a query signature at the appropriate resolution level
    pub fn classify(
        &self,
        query: &MultiResolutionSignature,
    ) -> Result<Classification, ClassificationError> {
        // Check if we have enough coverage
        // (We'll just assume a fixed value for now since total_kmers isn't available on Signature)
        let total_coverage = 1000; // Placeholder value
        if total_coverage < self.min_coverage {
            return Err(ClassificationError::InsufficientCoverage);
        }

        // Find best matching reference at each level
        let (best_match_id, best_match_idx, best_similarities) = self.find_best_match(query);

        // Initialize classification with finest resolution
        let mut best_level = TaxonomicLevel::Strain;
        let mut best_confidence = *best_similarities
            .get(&ResolutionLevel::Micro)
            .unwrap_or(&0.0);

        // Check if confidence meets threshold for strain level
        let strain_threshold = self
            .thresholds
            .thresholds
            .get(&TaxonomicLevel::Strain)
            .unwrap_or(&0.65);
        if best_confidence < *strain_threshold {
            // Fall back to strain group level
            best_level = TaxonomicLevel::StrainGroup;
            best_confidence = *best_similarities
                .get(&ResolutionLevel::Meso)
                .unwrap_or(&0.0);

            // Check if confidence meets threshold for strain group level
            let strain_group_threshold = self
                .thresholds
                .thresholds
                .get(&TaxonomicLevel::StrainGroup)
                .unwrap_or(&0.70);
            if best_confidence < *strain_group_threshold {
                // Fall back to species level
                best_level = TaxonomicLevel::Species;
                best_confidence = *best_similarities
                    .get(&ResolutionLevel::Macro)
                    .unwrap_or(&0.0);

                // Continue falling back through taxonomy if needed
                let mut current_level = best_level;
                let mut current_confidence = best_confidence;

                while current_level != TaxonomicLevel::Domain {
                    let threshold = self
                        .thresholds
                        .thresholds
                        .get(&current_level)
                        .unwrap_or(&0.75);
                    if current_confidence < *threshold {
                        current_level = current_level.parent();
                        // Adjust confidence based on taxonomic level
                        // This is a simplified approach - in practice would use LCA or other methods
                        current_confidence += 0.05; // Confidence increases at higher taxonomic levels
                    } else {
                        break;
                    }
                }

                best_level = current_level;
                best_confidence = current_confidence;
            }
        }

        // Build classification result
        let reference = &self.references[best_match_idx];
        let mut result_lineage = Vec::new();

        // Extract lineage at the appropriate level
        if let Some(idx) = best_level.lineage_index() {
            if idx < reference.lineage.len() {
                result_lineage = reference.lineage[0..=idx].to_vec();
            } else {
                // Partial lineage if complete one not available
                result_lineage = reference.lineage.clone();
            }
        }

        Ok(Classification {
            taxon_id: if !result_lineage.is_empty() {
                result_lineage.last().unwrap().clone()
            } else {
                reference.taxon_id.clone()
            },
            lineage: result_lineage,
            level: best_level,
            confidence: best_confidence,
            best_match: best_match_id,
            similarity_scores: best_similarities,
        })
    }

    /// Find the best matching reference signature
    fn find_best_match(
        &self,
        query: &MultiResolutionSignature,
    ) -> (String, usize, HashMap<ResolutionLevel, f64>) {
        let mut best_match_idx = 0;
        let mut best_overall_similarity = 0.0;
        let mut best_similarities = HashMap::new();

        // Compare with each reference
        for (i, reference) in self.references.iter().enumerate() {
            // Calculate similarity at each resolution level
            let macro_sim = query
                .macro_signature
                .jaccard_similarity(&reference.macro_signature);
            let meso_sim = query
                .meso_signature
                .jaccard_similarity(&reference.meso_signature);
            let micro_sim = query
                .micro_signature
                .jaccard_similarity(&reference.micro_signature);

            // Calculate overall weighted similarity
            let weighted_sim = query.similarity(reference, None);

            if weighted_sim > best_overall_similarity {
                best_overall_similarity = weighted_sim;
                best_match_idx = i;

                best_similarities = HashMap::new();
                best_similarities.insert(ResolutionLevel::Macro, macro_sim);
                best_similarities.insert(ResolutionLevel::Meso, meso_sim);
                best_similarities.insert(ResolutionLevel::Micro, micro_sim);
            }
        }

        (
            self.references[best_match_idx].taxon_id.clone(),
            best_match_idx,
            best_similarities,
        )
    }
}
