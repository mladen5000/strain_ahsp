use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};
use std::time::{Duration, SystemTime};

use bincode;
use log::{error, info, warn};
use rayon::prelude::*;
use reqwest::{blocking::Client, header};
use serde::{Deserialize, Serialize};
use sled::Db;
use thiserror::Error;

use crate::bio::signature::{MultiResolutionSignature, SignatureBuilder};

#[derive(Error, Debug)]
pub enum DatabaseError {
    #[error("IO error: {0}")]
    IoError(#[from] io::Error),

    #[error("HTTP error: {0}")]
    HttpError(#[from] reqwest::Error),

    #[error("NCBI API error: {0}")]
    NCBIApiError(String),

    #[error("Database error: {0}")]
    DatabaseError(#[from] sled::Error),

    #[error("Serialization error: {0}")]
    SerializationError(#[from] bincode::Error),

    #[error("Invalid taxonomy: {0}")]
    TaxonomyError(String),

    #[error("Signature error: {0}")]
    SignatureError(String),

    #[error("Not found: {0}")]
    NotFoundError(String),
}

/// NCBI genome metadata
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GenomeMetadata {
    /// NCBI accession number
    pub accession: String,

    /// Organism name
    pub organism: String,

    /// Taxonomy ID
    pub taxid: String,

    /// Assembly level (e.g., "Complete Genome", "Scaffold", etc.)
    pub assembly_level: String,

    /// Release date
    pub release_date: String,

    /// Size in base pairs
    pub size: usize,

    /// GC content
    pub gc_content: f64,

    /// Full taxonomic lineage
    pub lineage: Vec<(String, String)>, // (taxid, name) pairs
}

/// NCBI genome downloader
pub struct NCBIDownloader {
    /// HTTP client
    client: Client,

    /// Base URL for NCBI APIs
    base_url: String,

    /// API key for NCBI E-utilities (optional)
    api_key: Option<String>,

    /// Cache directory for downloaded genomes
    cache_dir: PathBuf,

    /// Cache expiration time in days
    cache_expiry_days: u64,
}

impl NCBIDownloader {
    /// Create a new NCBI downloader
    pub fn new(
        cache_dir: impl AsRef<Path>,
        api_key: Option<String>,
        cache_expiry_days: Option<u64>,
    ) -> Result<Self, DatabaseError> {
        // Create cache directory if it doesn't exist
        let cache_path = cache_dir.as_ref().to_path_buf();
        fs::create_dir_all(&cache_path)?;

        // Configure HTTP client with appropriate headers and timeouts
        let mut headers = header::HeaderMap::new();
        headers.insert(
            header::USER_AGENT,
            header::HeaderValue::from_static("ahsp/0.1 (https://github.com/your-repo/ahsp)"),
        );

        let client = Client::builder()
            .default_headers(headers)
            .timeout(Duration::from_secs(60))
            .build()?;

        Ok(NCBIDownloader {
            client,
            base_url: "https://eutils.ncbi.nlm.nih.gov/entrez/eutils".to_string(),
            api_key,
            cache_dir: cache_path,
            cache_expiry_days: cache_expiry_days.unwrap_or(30),
        })
    }

    /// Search for genomes matching a query
    pub fn search_genomes(
        &self,
        query: &str,
        max_results: usize,
    ) -> Result<Vec<GenomeMetadata>, DatabaseError> {
        // First, search for matching assembly IDs
        let search_url =
            format!(
            "{}/esearch.fcgi?db=assembly&term={}+AND+\"latest\"[filter]&retmax={}&retmode=json{}",
            self.base_url,
            urlencoding::encode(query),
            max_results,
            self.api_key.as_ref().map_or(String::new(), |k| format!("&api_key={}", k))
        );

        let search_response = self.client.get(&search_url).send()?;
        if !search_response.status().is_success() {
            return Err(DatabaseError::NCBIApiError(format!(
                "NCBI search failed with status: {}",
                search_response.status()
            )));
        }

        // Parse the search results to get assembly IDs
        let search_data: serde_json::Value = search_response.json()?;
        let id_list = match search_data["esearchresult"]["idlist"].as_array() {
            Some(list) if !list.is_empty() => list
                .iter()
                .filter_map(|id| id.as_str().map(|s| s.to_string()))
                .collect::<Vec<String>>(),
            _ => {
                return Err(DatabaseError::NCBIApiError(
                    "No assembly IDs found for the query".to_string(),
                ));
            }
        };

        // Now fetch details for each assembly ID
        let mut results = Vec::with_capacity(id_list.len());

        for id in id_list {
            // Fetch assembly summary
            let summary_url = format!(
                "{}/esummary.fcgi?db=assembly&id={}&retmode=json{}",
                self.base_url,
                id,
                self.api_key
                    .as_ref()
                    .map_or(String::new(), |k| format!("&api_key={}", k))
            );

            let summary_response = self.client.get(&summary_url).send()?;
            if !summary_response.status().is_success() {
                warn!(
                    "Failed to fetch summary for assembly ID {}: {}",
                    id,
                    summary_response.status()
                );
                continue;
            }

            let summary_data: serde_json::Value = summary_response.json()?;
            let result = &summary_data["result"][&id];

            // Extract metadata
            let accession = result["assemblyaccession"]
                .as_str()
                .unwrap_or("")
                .to_string();

            let organism = result["speciesname"].as_str().unwrap_or("").to_string();

            let taxid = result["taxid"].as_str().unwrap_or("").to_string();

            let assembly_level = result["assemblylevel"].as_str().unwrap_or("").to_string();

            let release_date = result["submissiondate"].as_str().unwrap_or("").to_string();

            let size = result["totallength"]
                .as_str()
                .unwrap_or("0")
                .parse::<usize>()
                .unwrap_or(0);

            let gc_content = result["genomegcpercent"].as_f64().unwrap_or(0.0);

            // Fetch taxonomy lineage
            let lineage = self.fetch_taxonomy_lineage(&taxid)?;

            results.push(GenomeMetadata {
                accession,
                organism,
                taxid,
                assembly_level,
                release_date,
                size,
                gc_content,
                lineage,
            });
        }

        Ok(results)
    }

    /// Fetch taxonomic lineage for a taxonomy ID
    fn fetch_taxonomy_lineage(&self, taxid: &str) -> Result<Vec<(String, String)>, DatabaseError> {
        // Fetch full taxonomy information
        let elink_url = format!(
            "{}/efetch.fcgi?db=taxonomy&id={}&retmode=xml{}",
            self.base_url,
            taxid,
            self.api_key
                .as_ref()
                .map_or(String::new(), |k| format!("&api_key={}", k))
        );

        let response = self.client.get(&elink_url).send()?;
        if !response.status().is_success() {
            return Err(DatabaseError::NCBIApiError(format!(
                "Taxonomy fetch failed with status: {}",
                response.status()
            )));
        }

        // Parse XML response to extract lineage
        let xml_text = response.text()?;
        let mut lineage = Vec::new();

        // Use quick-xml to parse the taxonomy lineage
        let mut reader = quick_xml::Reader::from_str(&xml_text);
        let mut buf = Vec::new();

        let mut current_name = String::new();
        let mut current_id = String::new();
        let mut in_lineage = false;
        let mut in_taxid = false;
        let mut in_name = false;

        loop {
            match reader.read_event(&mut buf) {
                Ok(quick_xml::events::Event::Start(ref e)) => match e.name() {
                    b"Lineage" => in_lineage = true,
                    b"TaxId" => in_taxid = true,
                    b"ScientificName" => in_name = true,
                    _ => {}
                },
                Ok(quick_xml::events::Event::End(ref e)) => match e.name() {
                    b"Lineage" => in_lineage = false,
                    b"TaxId" => in_taxid = false,
                    b"ScientificName" => in_name = false,
                    _ => {}
                },
                Ok(quick_xml::events::Event::Text(ref e)) => {
                    if in_taxid {
                        current_id = e.unescape_and_decode(&reader).unwrap_or_default();
                    } else if in_name {
                        current_name = e.unescape_and_decode(&reader).unwrap_or_default();
                    } else if in_lineage {
                        // Process lineage text - typically a semicolon-separated list
                        let lineage_text = e.unescape_and_decode(&reader).unwrap_or_default();
                        for name in lineage_text.split(';') {
                            let name = name.trim();
                            if !name.is_empty() {
                                // We don't have taxids for these lineage entries,
                                // so we use empty strings as placeholders
                                lineage.push(("".to_string(), name.to_string()));
                            }
                        }
                    }
                }
                Ok(quick_xml::events::Event::Eof) => break,
                Err(e) => {
                    return Err(DatabaseError::TaxonomyError(format!(
                        "XML parsing error: {}",
                        e
                    )));
                }
                _ => {}
            }
            buf.clear();
        }

        // Add the current taxon at the end of the lineage
        if !current_id.is_empty() && !current_name.is_empty() {
            lineage.push((current_id, current_name));
        }

        Ok(lineage)
    }

    /// Download a genome FASTA file by accession
    pub fn download_genome(&self, accession: &str) -> Result<PathBuf, DatabaseError> {
        // Check if genome is already cached
        let cache_file = self.cache_dir.join(format!("{}.fna.gz", accession));

        // Check if cache file exists and is not expired
        if cache_file.exists() {
            if let Ok(metadata) = fs::metadata(&cache_file) {
                if let Ok(modified) = metadata.modified() {
                    let now = SystemTime::now();
                    if now
                        .duration_since(modified)
                        .map(|d| d.as_secs() < self.cache_expiry_days * 86400)
                        .unwrap_or(false)
                    {
                        // Cache is still valid
                        return Ok(cache_file);
                    }
                }
            }
        }

        // Cache doesn't exist or is expired, download the genome
        info!("Downloading genome for accession: {}", accession);

        // First, we need to find the FTP path for this assembly
        let summary_url = format!(
            "{}/esummary.fcgi?db=assembly&id={}&retmode=json{}",
            self.base_url,
            accession,
            self.api_key
                .as_ref()
                .map_or(String::new(), |k| format!("&api_key={}", k))
        );

        let summary_response = self.client.get(&summary_url).send()?;
        if !summary_response.status().is_success() {
            return Err(DatabaseError::NCBIApiError(format!(
                "Assembly summary fetch failed with status: {}",
                summary_response.status()
            )));
        }

        let summary_data: serde_json::Value = summary_response.json()?;
        let ftp_path = summary_data["result"][accession]["ftppath_genbank"]
            .as_str()
            .or_else(|| summary_data["result"][accession]["ftppath_refseq"].as_str())
            .ok_or_else(|| {
                DatabaseError::NCBIApiError("No FTP path found for assembly".to_string())
            })?;

        // Construct the FASTA download URL
        let filename = format!(
            "{}_{}.fna.gz",
            accession,
            summary_data["result"][accession]["assemblyname"]
                .as_str()
                .unwrap_or("genomic")
                .replace(" ", "_")
        );

        let download_url = format!("{}/{}", ftp_path, filename);

        // Download the file
        let response = self.client.get(&download_url).send()?;
        if !response.status().is_success() {
            return Err(DatabaseError::NCBIApiError(format!(
                "Genome download failed with status: {}",
                response.status()
            )));
        }

        // Save to cache
        let mut file = File::create(&cache_file)?;
        let content = response.bytes()?;
        file.write_all(&content)?;

        Ok(cache_file)
    }
}

/// Signature database
pub struct SignatureDatabase {
    /// Underlying key-value store
    db: Db,

    /// Index of taxonomy IDs to accessions
    taxonomy_index: HashMap<String, Vec<String>>,

    /// Index of lineage terms to accessions
    lineage_index: HashMap<String, Vec<String>>,
}

impl SignatureDatabase {
    /// Open or create a signature database
    pub fn open(path: impl AsRef<Path>) -> Result<Self, DatabaseError> {
        // Open or create sled database
        let db = sled::open(path)?;

        // Load indices
        let taxonomy_index = match db.get("taxonomy_index")? {
            Some(data) => bincode::deserialize(&data)?,
            None => HashMap::new(),
        };

        let lineage_index = match db.get("lineage_index")? {
            Some(data) => bincode::deserialize(&data)?,
            None => HashMap::new(),
        };

        Ok(SignatureDatabase {
            db,
            taxonomy_index,
            lineage_index,
        })
    }

    /// Add a signature to the database
    pub fn add_signature(
        &mut self,
        signature: &MultiResolutionSignature,
    ) -> Result<(), DatabaseError> {
        // Serialize the signature
        let signature_data = bincode::serialize(signature)?;

        // Store in database
        self.db.insert(&signature.taxon_id, signature_data)?;

        // Update indices
        self.update_indices(signature)?;

        Ok(())
    }

    /// Update search indices for a signature
    fn update_indices(
        &mut self,
        signature: &MultiResolutionSignature,
    ) -> Result<(), DatabaseError> {
        // Update taxonomy index
        self.taxonomy_index
            .entry(signature.taxon_id.clone())
            .or_insert_with(Vec::new)
            .push(signature.taxon_id.clone());

        // Update lineage index
        for taxon in &signature.lineage {
            self.lineage_index
                .entry(taxon.clone())
                .or_insert_with(Vec::new)
                .push(signature.taxon_id.clone());
        }

        // Save indices
        let taxonomy_data = bincode::serialize(&self.taxonomy_index)?;
        let lineage_data = bincode::serialize(&self.lineage_index)?;

        self.db.insert("taxonomy_index", taxonomy_data)?;
        self.db.insert("lineage_index", lineage_data)?;

        Ok(())
    }

    /// Get a signature by ID
    pub fn get_signature(&self, id: &str) -> Result<MultiResolutionSignature, DatabaseError> {
        match self.db.get(id)? {
            Some(data) => Ok(bincode::deserialize(&data)?),
            None => Err(DatabaseError::NotFoundError(format!(
                "Signature not found: {}",
                id
            ))),
        }
    }

    /// Search for signatures matching a taxonomy term
    pub fn search_by_taxonomy(
        &self,
        term: &str,
    ) -> Result<Vec<MultiResolutionSignature>, DatabaseError> {
        let mut results = Vec::new();

        // Search taxonomy index
        if let Some(ids) = self.taxonomy_index.get(term) {
            for id in ids {
                if let Ok(signature) = self.get_signature(id) {
                    results.push(signature);
                }
            }
        }

        // Search lineage index
        if let Some(ids) = self.lineage_index.get(term) {
            for id in ids {
                if let Ok(signature) = self.get_signature(id) {
                    if !results.iter().any(|s| s.taxon_id == *id) {
                        results.push(signature);
                    }
                }
            }
        }

        Ok(results)
    }

    /// Get all signatures
    pub fn get_all_signatures(&self) -> Result<Vec<MultiResolutionSignature>, DatabaseError> {
        let mut results = Vec::new();

        for item in self.db.iter() {
            let (key, value) = item?;
            let key_str = std::str::from_utf8(&key)
                .map_err(|_| DatabaseError::DatabaseError("Invalid UTF-8 in key".into()))?;

            // Skip index entries
            if key_str != "taxonomy_index" && key_str != "lineage_index" {
                let signature: MultiResolutionSignature = bincode::deserialize(&value)?;
                results.push(signature);
            }
        }

        Ok(results)
    }
}

/// Database manager for coordinating NCBI downloads and signature generation
pub struct DatabaseManager {
    /// Signature database
    pub database: SignatureDatabase,

    /// NCBI downloader
    pub downloader: NCBIDownloader,

    /// Signature builder
    pub builder: SignatureBuilder,
}

impl DatabaseManager {
    /// Create a new database manager
    pub fn new(
        db_path: impl AsRef<Path>,
        cache_dir: impl AsRef<Path>,
        macro_k: usize,
        meso_k: usize,
        sketch_size: usize,
        threads: usize,
        api_key: Option<String>,
    ) -> Result<Self, DatabaseError> {
        let database = SignatureDatabase::open(db_path)?;
        let downloader = NCBIDownloader::new(cache_dir, api_key, None)?;

        let builder = SignatureBuilder::new(macro_k, meso_k, sketch_size, threads)
            .map_err(|e| DatabaseError::SignatureError(format!("{}", e)))?;

        Ok(DatabaseManager {
            database,
            downloader,
            builder,
        })
    }

    /// Search for and download reference genomes from NCBI
    pub fn download_references(
        &self,
        query: &str,
        max_results: usize,
    ) -> Result<Vec<(GenomeMetadata, PathBuf)>, DatabaseError> {
        // Search for matching genomes
        let genomes = self.downloader.search_genomes(query, max_results)?;

        // Download each genome
        let results = genomes
            .into_par_iter()
            .filter_map(
                |genome| match self.downloader.download_genome(&genome.accession) {
                    Ok(path) => Some((genome, path)),
                    Err(e) => {
                        error!("Failed to download genome {}: {}", genome.accession, e);
                        None
                    }
                },
            )
            .collect();

        Ok(results)
    }

    /// Process downloaded genomes into signatures and add to database
    pub fn process_references(
        &mut self,
        references: Vec<(GenomeMetadata, PathBuf)>,
    ) -> Result<Vec<String>, DatabaseError> {
        // Convert GenomeMetadata to format needed by SignatureBuilder
        let files = references
            .iter()
            .map(|(metadata, path)| {
                // Convert lineage to string vector for signature
                let lineage = metadata
                    .lineage
                    .iter()
                    .map(|(_, name)| name.clone())
                    .collect();

                (path.clone(), metadata.accession.clone(), lineage)
            })
            .collect::<Vec<_>>();

        // Build signatures in batch
        let signatures = self
            .builder
            .build_batch(files)
            .map_err(|e| DatabaseError::SignatureError(format!("{}", e)))?;

        // Add signatures to database
        let mut added_ids = Vec::with_capacity(signatures.len());
        for signature in signatures {
            self.database.add_signature(&signature)?;
            added_ids.push(signature.taxon_id.clone());
        }

        Ok(added_ids)
    }

    /// Search and add reference genomes by query
    pub fn search_and_add_references(
        &mut self,
        query: &str,
        max_results: usize,
    ) -> Result<Vec<String>, DatabaseError> {
        // Download references
        let references = self.download_references(query, max_results)?;

        // Process references
        self.process_references(references)
    }

    /// Check if the database is empty
    pub fn is_empty(&self) -> Result<bool, DatabaseError> {
        let count = self.database.get_all_signatures()?.len();
        Ok(count == 0)
    }
}
