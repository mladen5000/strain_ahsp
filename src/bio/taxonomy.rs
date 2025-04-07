//! Taxonomic classification utilities.
//!
//! This module provides structures and functions for working with
//! taxonomic classifications and lineages.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Taxonomic classification levels.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum TaxonomicLevel {
    Domain,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species,
    Strain,
}

impl TaxonomicLevel {
    /// Returns a string representation of the taxonomic level.
    pub fn as_str(&self) -> &'static str {
        match self {
            TaxonomicLevel::Domain => "domain",
            TaxonomicLevel::Phylum => "phylum",
            TaxonomicLevel::Class => "class",
            TaxonomicLevel::Order => "order",
            TaxonomicLevel::Family => "family",
            TaxonomicLevel::Genus => "genus",
            TaxonomicLevel::Species => "species",
            TaxonomicLevel::Strain => "strain",
        }
    }
    
    /// Returns the hierarchical depth of this level.
    pub fn depth(&self) -> usize {
        match self {
            TaxonomicLevel::Domain => 1,
            TaxonomicLevel::Phylum => 2,
            TaxonomicLevel::Class => 3,
            TaxonomicLevel::Order => 4,
            TaxonomicLevel::Family => 5,
            TaxonomicLevel::Genus => 6,
            TaxonomicLevel::Species => 7,
            TaxonomicLevel::Strain => 8,
        }
    }
    
    /// Returns all taxonomic levels in hierarchical order.
    pub fn all_levels() -> Vec<TaxonomicLevel> {
        vec![
            TaxonomicLevel::Domain,
            TaxonomicLevel::Phylum,
            TaxonomicLevel::Class,
            TaxonomicLevel::Order,
            TaxonomicLevel::Family,
            TaxonomicLevel::Genus,
            TaxonomicLevel::Species,
            TaxonomicLevel::Strain,
        ]
    }
}

/// Represents a complete taxonomic lineage from domain to strain.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct TaxonomicLineage {
    // Maps taxonomic levels to their taxon names
    levels: HashMap<TaxonomicLevel, String>,
    // Optional taxonomic ID (e.g., NCBI taxon ID)
    tax_id: Option<String>,
}

impl TaxonomicLineage {
    /// Creates a new empty taxonomic lineage.
    pub fn new() -> Self {
        TaxonomicLineage {
            levels: HashMap::new(),
            tax_id: None,
        }
    }
    
    /// Creates a new taxonomic lineage with the given taxonomy ID.
    pub fn with_tax_id(tax_id: String) -> Self {
        TaxonomicLineage {
            levels: HashMap::new(),
            tax_id: Some(tax_id),
        }
    }
    
    /// Sets a taxonomic level with its taxon name.
    pub fn set_level(&mut self, level: TaxonomicLevel, name: String) {
        self.levels.insert(level, name);
    }
    
    /// Gets the taxon name at a specific taxonomic level.
    pub fn get_level(&self, level: TaxonomicLevel) -> Option<&String> {
        self.levels.get(&level)
    }
    
    /// Returns the most specific (deepest) taxonomic level with a defined name.
    pub fn most_specific_level(&self) -> Option<TaxonomicLevel> {
        TaxonomicLevel::all_levels()
            .into_iter()
            .rev() // Start from strain and go up
            .find(|level| self.levels.contains_key(level))
    }
    
    /// Returns a vector of all defined taxonomic levels and their names.
    pub fn to_vec(&self) -> Vec<(TaxonomicLevel, &String)> {
        let mut result: Vec<(TaxonomicLevel, &String)> = self.levels
            .iter()
            .map(|(level, name)| (*level, name))
            .collect();
        
        // Sort by taxonomic depth
        result.sort_by_key(|(level, _)| level.depth());
        
        result
    }
    
    /// Returns the taxonomy ID, if set.
    pub fn tax_id(&self) -> Option<&String> {
        self.tax_id.as_ref()
    }
    
    /// Sets the taxonomy ID.
    pub fn set_tax_id(&mut self, tax_id: String) {
        self.tax_id = Some(tax_id);
    }
    
    /// Returns a string representation of the lineage, from domain to the most specific level.
    pub fn to_string(&self) -> String {
        let parts: Vec<String> = self.to_vec()
            .into_iter()
            .map(|(_, name)| name.clone())
            .collect();
        
        parts.join("; ")
    }
}

/// Parses a taxonomic lineage from a string in standard format.
/// 
/// # Arguments
/// 
/// * `lineage_str` - A string like "Bacteria; Proteobacteria; Gammaproteobacteria; ..."
/// 
/// # Returns
/// 
/// A TaxonomicLineage with levels parsed from the string
pub fn parse_lineage(lineage_str: &str) -> TaxonomicLineage {
    let parts: Vec<&str> = lineage_str.split(';').map(str::trim).collect();
    let mut lineage = TaxonomicLineage::new();
    
    let levels = TaxonomicLevel::all_levels();
    
    for (i, part) in parts.iter().enumerate() {
        if i < levels.len() && !part.is_empty() {
            lineage.set_level(levels[i], part.to_string());
        }
    }
    
    lineage
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_taxonomic_level_as_str() {
        assert_eq!(TaxonomicLevel::Domain.as_str(), "domain");
        assert_eq!(TaxonomicLevel::Species.as_str(), "species");
    }
    
    #[test]
    fn test_taxonomic_level_depth() {
        assert_eq!(TaxonomicLevel::Domain.depth(), 1);
        assert_eq!(TaxonomicLevel::Species.depth(), 7);
    }
    
    #[test]
    fn test_taxonomic_lineage_basics() {
        let mut lineage = TaxonomicLineage::new();
        
        // Initially empty
        assert!(lineage.get_level(TaxonomicLevel::Domain).is_none());
        assert_eq!(lineage.to_string(), "");
        
        // Add some levels
        lineage.set_level(TaxonomicLevel::Domain, "Bacteria".to_string());
        lineage.set_level(TaxonomicLevel::Phylum, "Proteobacteria".to_string());
        
        // Check levels
        assert_eq!(lineage.get_level(TaxonomicLevel::Domain).unwrap(), "Bacteria");
        assert_eq!(lineage.to_string(), "Bacteria; Proteobacteria");
        
        // Most specific level
        assert_eq!(lineage.most_specific_level(), Some(TaxonomicLevel::Phylum));
        
        // Check vector representation
        let vec_rep = lineage.to_vec();
        assert_eq!(vec_rep.len(), 2);
        assert_eq!(vec_rep[0].0, TaxonomicLevel::Domain);
        assert_eq!(vec_rep[1].0, TaxonomicLevel::Phylum);
    }
    
    #[test]
    fn test_parse_lineage() {
        let lineage_str = "Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales; Enterobacteriaceae; Escherichia; Escherichia coli";
        let lineage = parse_lineage(lineage_str);
        
        assert_eq!(lineage.get_level(TaxonomicLevel::Domain).unwrap(), "Bacteria");
        assert_eq!(lineage.get_level(TaxonomicLevel::Genus).unwrap(), "Escherichia");
        assert_eq!(lineage.get_level(TaxonomicLevel::Species).unwrap(), "Escherichia coli");
        assert_eq!(lineage.most_specific_level(), Some(TaxonomicLevel::Species));
    }
    
    #[test]
    fn test_lineage_with_tax_id() {
        let mut lineage = TaxonomicLineage::with_tax_id("562".to_string()); // E. coli NCBI taxon ID
        
        assert_eq!(lineage.tax_id().unwrap(), "562");
        
        lineage.set_tax_id("1280".to_string()); // Changed to another ID
        assert_eq!(lineage.tax_id().unwrap(), "1280");
    }
}