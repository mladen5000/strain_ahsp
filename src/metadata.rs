//! Metadata handling module.
//!
//! This module provides structures and functions for working with sample metadata,
//! including experimental design and sample information.

use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;

#[derive(Debug, Serialize, Deserialize)]
pub struct Metadata {
    pub sample_info: HashMap<String, SampleInfo>,
    pub condition_map: HashMap<String, String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct SampleInfo {
    pub condition: String,
    pub replicate: u32,
    // Add other metadata fields as needed
}

pub fn load_metadata(path: &str) -> Result<Metadata> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let metadata = csv::Reader::from_reader(reader)
        .into_deserialize()
        .collect::<Result<Vec<(String, SampleInfo)>, _>>()?;
    Ok(Metadata {
        sample_info: metadata.into_iter().collect(),
        condition_map: HashMap::new(),
    })
}
