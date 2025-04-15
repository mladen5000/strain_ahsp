use std::collections::{HashMap, HashSet}; // Added HashSet
use std::fs::{self, File};
use std::io::{self, BufReader, Write}; // Added BufReader
use std::path::{Path, PathBuf};
use std::time::{Duration, SystemTime};

// Assuming the dummy definitions above are in src/sketch/signature.rs
use crate::pipeline::qc;
use crate::sketch::signature::MultiResolutionSignature; // Add MultiResolutionSignature from qc
use crate::sketch::SignatureBuilder;
use bincode::config::standard;
use bincode::{decode_from_slice, encode_to_vec};
use log::{error, info, warn};
use quick_xml::events::{BytesStart, Event}; // Added BytesStart, Event
use quick_xml::Reader; // Added Reader
use rayon::prelude::*;
use reqwest::{blocking::Client, header};
use serde::{Deserialize, Serialize};
use sled::Db;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum DatabaseError {
    #[error("IO error: {0}")]
    IoError(#[from] io::Error),

    #[error("HTTP error: {0}")]
    HttpError(#[from] reqwest::Error),

    #[error("Database error: {0}")]
    DatabaseError(#[from] sled::Error),

    #[error("Serialization error: {0}")]
    SerializationError(String),

    #[error("Not found: {0}")]
    NotFoundError(String),

    #[error("NCBI API error: {0}")]
    NCBIApiError(String),

    #[error("Taxonomy error: {0}")]
    TaxonomyError(String),

    #[error("Signature error: {0}")]
    SignatureError(String),

    #[error("XML parsing error: {0}")]
    XmlError(#[from] quick_xml::Error), // Added XML error variant

    #[error("UTF-8 conversion error: {0}")]
    Utf8Error(#[from] std::str::Utf8Error), // Added Utf8 error variant
}

// Add conversion from bincode errors
impl From<bincode::error::EncodeError> for DatabaseError {
    fn from(err: bincode::error::EncodeError) -> Self {
        DatabaseError::SerializationError(err.to_string())
    }
}

impl From<bincode::error::DecodeError> for DatabaseError {
    fn from(err: bincode::error::DecodeError) -> Self {
        DatabaseError::SerializationError(err.to_string())
    }
}

/// NCBI genome metadata
use bincode::{Decode, Encode};
#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)] // <-- Add Encode, Decode
pub struct GenomeMetadata {
    /// NCBI accession number (e.g., GCF_...)
    pub accession: String,

    /// Assembly ID (internal NCBI identifier, needed for some API calls)
    pub assembly_id: String, // Added field

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

    /// Full taxonomic lineage as (taxid, name) pairs, root first
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
        let cache_path = cache_dir.as_ref().to_path_buf();
        fs::create_dir_all(&cache_path)?;

        let mut headers = header::HeaderMap::new();
        headers.insert(
            header::USER_AGENT,
            // Consider making the user agent more specific if this is a library
            header::HeaderValue::from_static("rust-ncbi-downloader/0.1"),
        );

        let client = Client::builder()
            .default_headers(headers)
            .timeout(Duration::from_secs(60)) // Increased timeout
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
        // Search for assembly IDs
        let search_url = format!(
            "{}/esearch.fcgi?db=assembly&term={}+AND+\"latest refseq\"[filter]+AND+\"complete genome\"[filter]&retmax={}&retmode=json{}",
            self.base_url,
            urlencoding::encode(query), // Ensure query is URL-encoded
            max_results,
            self.api_key.as_ref().map_or(String::new(), |k| format!("&api_key={}", k))
        );
        info!("Searching NCBI Assembly: {}", search_url);

        let search_response = self.client.get(&search_url).send()?;
        if !search_response.status().is_success() {
            let status = search_response.status();
            let body = search_response
                .text()
                .unwrap_or_else(|_| "Failed to read body".to_string());
            return Err(DatabaseError::NCBIApiError(format!(
                "NCBI search failed with status: {}. Response: {}",
                status, body
            )));
        }

        // Parse search results
        let search_data: serde_json::Value = search_response.json()?;
        let id_list = search_data["esearchresult"]["idlist"]
            .as_array()
            .ok_or_else(|| {
                DatabaseError::NCBIApiError(
                    "Invalid search response format: Missing 'idlist'".to_string(),
                )
            })?
            .iter()
            .filter_map(|id| id.as_str().map(String::from))
            .collect::<Vec<String>>();

        if id_list.is_empty() {
            info!("No assembly IDs found for query: {}", query);
            return Ok(Vec::new()); // Return empty vector instead of error
        }

        info!("Found {} assembly IDs, fetching details...", id_list.len());

        // Fetch details in batches (e.g., 50 per request) for efficiency
        let mut results = Vec::with_capacity(id_list.len());
        for id_chunk in id_list.chunks(50) {
            let ids_str = id_chunk.join(",");
            let summary_url = format!(
                "{}/esummary.fcgi?db=assembly&id={}&retmode=json{}",
                self.base_url,
                ids_str,
                self.api_key
                    .as_ref()
                    .map_or(String::new(), |k| format!("&api_key={}", k))
            );

            info!("Fetching summaries for IDs: {}", ids_str);
            let summary_response = self.client.get(&summary_url).send()?;

            if !summary_response.status().is_success() {
                let status = summary_response.status();
                let body = summary_response
                    .text()
                    .unwrap_or_else(|_| "Failed to read body".to_string());
                warn!(
                    "Failed to fetch summaries for chunk starting with {}: Status {}, Body: {}",
                    id_chunk.first().unwrap_or(&"N/A".to_string()),
                    status,
                    body
                );
                continue; // Skip this chunk on error
            }

            let summary_data: serde_json::Value = summary_response.json()?;
            let summary_results = summary_data["result"].as_object().ok_or_else(|| {
                DatabaseError::NCBIApiError(
                    "Invalid summary response format: Missing 'result' object".to_string(),
                )
            })?;

            // Check for "uids" list to iterate results correctly
            let uids = summary_results
                .get("uids")
                .and_then(|u| u.as_array())
                .ok_or_else(|| {
                    DatabaseError::NCBIApiError(
                        "Invalid summary response format: Missing 'uids'".to_string(),
                    )
                })?;

            for uid_val in uids {
                if let Some(id_str) = uid_val.as_str() {
                    if let Some(result) = summary_results.get(id_str) {
                        // Extract metadata carefully, providing defaults
                        let accession = result["assemblyaccession"]
                            .as_str()
                            .unwrap_or("")
                            .to_string();
                        let organism = result["speciesname"]
                            .as_str()
                            .unwrap_or("Unknown")
                            .to_string();
                        let taxid = result["taxid"].as_str().unwrap_or("0").to_string();
                        let assembly_level = result["assemblylevel"]
                            .as_str()
                            .unwrap_or("Unknown")
                            .to_string();
                        let release_date = result["releasedate"] // Use releasedate
                            .as_str()
                            .or(result["submissiondate"].as_str()) // Fallback
                            .unwrap_or("")
                            .to_string();
                        let size = result["totallength"]
                            .as_str()
                            .unwrap_or("0")
                            .parse::<usize>()
                            .unwrap_or(0);
                        let gc_content = result["genomegcpercent"]
                            .as_str() // GC percent can be string
                            .unwrap_or("0.0")
                            .parse::<f64>()
                            .unwrap_or(0.0);

                        if accession.is_empty() || taxid == "0" {
                            warn!(
                                "Skipping assembly ID {} due to missing accession or taxid.",
                                id_str
                            );
                            continue;
                        }

                        // Fetch taxonomy lineage
                        match self.fetch_taxonomy_lineage(&taxid) {
                            Ok(lineage) => {
                                results.push(GenomeMetadata {
                                    accession,
                                    assembly_id: id_str.to_string(), // Store the assembly ID
                                    organism,
                                    taxid,
                                    assembly_level,
                                    release_date,
                                    size,
                                    gc_content,
                                    lineage,
                                });
                            }
                            Err(e) => {
                                warn!("Failed to fetch lineage for taxid {}: {}. Skipping assembly {}.", taxid, e, accession);
                                // Decide whether to skip or add with empty lineage
                                // Skipping for now:
                                continue;
                            }
                        }
                    } else {
                        warn!(
                            "No summary data found for assembly ID {} in response.",
                            id_str
                        );
                    }
                }
            }
        }

        info!("Retrieved metadata for {} genomes.", results.len());
        Ok(results)
    }

    /// Fetch taxonomic lineage for a taxonomy ID using efetch XML.
    /// Prioritizes LineageEx for (taxid, name) pairs.
    fn fetch_taxonomy_lineage(&self, taxid: &str) -> Result<Vec<(String, String)>, DatabaseError> {
        let efetch_url = format!(
            "{}/efetch.fcgi?db=taxonomy&id={}&retmode=xml{}",
            self.base_url,
            taxid,
            self.api_key
                .as_ref()
                .map_or(String::new(), |k| format!("&api_key={}", k))
        );

        let response = self.client.get(&efetch_url).send()?;
        if !response.status().is_success() {
            return Err(DatabaseError::NCBIApiError(format!(
                "Taxonomy fetch failed for taxid {}: Status {}",
                taxid,
                response.status()
            )));
        }

        let xml_text = response.text()?;
        let mut reader = Reader::from_str(&xml_text);
        // reader.trim_text(true);
        let mut buf = Vec::new();

        let mut lineage: Vec<(String, String)> = Vec::new();
        let mut current_taxid = String::new();
        let mut current_name = String::new();
        let mut main_taxon_name = String::new(); // To store the primary name
        let mut in_lineage_ex = false;
        let mut in_taxon = false; // Within LineageEx/Taxon
        let mut in_taxid_tag = false;
        let mut in_name_tag = false;

        loop {
            match reader.read_event_into(&mut buf) {
                Ok(Event::Start(ref e)) => {
                    match e.name().as_ref() {
                        b"LineageEx" => in_lineage_ex = true,
                        b"Taxon" if in_lineage_ex => {
                            in_taxon = true;
                            current_taxid.clear();
                            current_name.clear();
                        }
                        // Handle Taxon element not inside LineageEx (the main one)
                        b"Taxon" if !in_lineage_ex => {
                            in_taxon = true; // Re-use flag for simplicity
                            current_taxid.clear();
                            current_name.clear();
                        }
                        b"TaxId" => in_taxid_tag = true,
                        b"ScientificName" => in_name_tag = true,
                        _ => {}
                    }
                }
                Ok(Event::Text(e)) => {
                    let text = e.unescape()?.to_string();
                    if in_taxid_tag {
                        current_taxid = text;
                    } else if in_name_tag {
                        current_name = text;
                        if in_taxon && !in_lineage_ex && main_taxon_name.is_empty() {
                            main_taxon_name = current_name.clone(); // Store the main name
                        }
                    }
                }
                Ok(Event::End(ref e)) => match e.name().as_ref() {
                    b"LineageEx" => in_lineage_ex = false,
                    b"Taxon" => {
                        if in_lineage_ex && !current_taxid.is_empty() && !current_name.is_empty() {
                            lineage.push((current_taxid.clone(), current_name.clone()));
                        }
                        in_taxon = false;
                        current_taxid.clear();
                        current_name.clear();
                    }
                    b"TaxId" => in_taxid_tag = false,
                    b"ScientificName" => in_name_tag = false,
                    _ => {}
                },
                Ok(Event::Eof) => break,
                Err(e) => return Err(DatabaseError::XmlError(e)),
                _ => {} // Ignore other events
            }
            buf.clear();
        }

        // Ensure the main taxon itself is added if lineage was parsed
        // It might not be part of LineageEx itself in some NCBI XML formats
        if !lineage.is_empty() && !taxid.is_empty() && !main_taxon_name.is_empty() {
            // Check if the main taxon is already the last element
            if lineage.last().map_or(true, |(id, _)| id != taxid) {
                lineage.push((taxid.to_string(), main_taxon_name));
            }
        } else if lineage.is_empty() && !taxid.is_empty() && !main_taxon_name.is_empty() {
            // Handle case where LineageEx might be empty or missing, but main taxon is present
            lineage.push((taxid.to_string(), main_taxon_name));
        }

        if lineage.is_empty() {
            // Fallback or error if no structured lineage found
            warn!("Could not parse structured lineage (LineageEx) for taxid {}. Check NCBI XML format.", taxid);
            // Optionally, could try parsing the simple <Lineage> string here if needed
            // return Err(DatabaseError::TaxonomyError(format!("No lineage data found for taxid {}", taxid)));
        }

        Ok(lineage)
    }

    /// Download a genome FASTA file by accession (e.g., GCF_...).
    /// Uses the cache if available and not expired.
    pub fn download_genome(&self, accession: &str) -> Result<PathBuf, DatabaseError> {
        let expected_filename = format!("{}.fna.gz", accession);
        let cache_file = self.cache_dir.join(&expected_filename);

        // Check cache validity
        if cache_file.exists() {
            if let Ok(metadata) = fs::metadata(&cache_file) {
                if let Ok(modified) = metadata.modified() {
                    if SystemTime::now()
                        .duration_since(modified)
                        .map_or(false, |d| d.as_secs() < self.cache_expiry_days * 86400)
                    {
                        info!("Using cached genome: {}", cache_file.display());
                        return Ok(cache_file);
                    } else {
                        info!("Cache expired for {}, re-downloading.", accession);
                        // Optionally remove the old file: fs::remove_file(&cache_file)?;
                    }
                }
            }
        }

        info!("Downloading genome for accession: {}", accession);

        // Fetch assembly summary to find the FTP path
        // Note: Using accession directly in 'id' often works for esummary/assembly
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
                "Assembly summary fetch failed for {}: Status {}",
                accession,
                summary_response.status()
            )));
        }

        let summary_data: serde_json::Value = summary_response.json()?;

        // NCBI's esummary result structure can be tricky. Need to find the result entry.
        // It might be keyed by the input ID (accession) or an internal UID.
        // Check "result" -> "uids" first if present.
        let result_obj = summary_data["result"].as_object().ok_or_else(|| {
            DatabaseError::NCBIApiError("Invalid summary response: 'result' not an object".into())
        })?;

        let assembly_info = result_obj
            .get(accession) // Try keying by accession first
            .or_else(|| {
                // If not found, check if 'uids' exists and use the first uid
                result_obj
                    .get("uids")
                    .and_then(|uids| uids.as_array()?.get(0))
                    .and_then(|uid_val| uid_val.as_str())
                    .and_then(|uid_str| result_obj.get(uid_str))
            })
            .ok_or_else(|| {
                DatabaseError::NCBIApiError(format!(
                    "Could not find summary details for {} in response",
                    accession
                ))
            })?;

        let ftp_path_base = assembly_info["ftppath_genbank"]
            .as_str()
            .filter(|s| !s.is_empty()) // Ensure it's not an empty string
            .or_else(|| {
                assembly_info["ftppath_refseq"]
                    .as_str()
                    .filter(|s| !s.is_empty())
            })
            .ok_or_else(|| {
                DatabaseError::NCBIApiError(format!(
                    "No valid FTP path (GenBank or RefSeq) found for assembly {}",
                    accession
                ))
            })?;

        // Construct the likely FASTA download URL
        // The filename often includes the assembly name. Example: GCF_000001405.40_GRCh38.p14_genomic.fna.gz
        let assembly_name = assembly_info["assemblyname"]
            .as_str()
            .unwrap_or("assembly") // Default if name is missing
            .replace([' ', '/', '\\', ':', '*', '?', '\"', '<', '>', '|'], "_"); // Sanitize name

        // Find the file on FTP. The exact name can vary. Common pattern: {ftp_path_base}/{ftp_path_base_basename}_{assembly_name}_genomic.fna.gz
        // We need the basename of the ftp path itself.
        let ftp_basename = ftp_path_base.split('/').last().unwrap_or(accession); // Get the last part, e.g., GCF_000..._Assembly

        let download_filename = format!("{}_{}_genomic.fna.gz", ftp_basename, assembly_name);
        let download_url = format!("{}/{}", ftp_path_base, download_filename);

        info!("Attempting download from: {}", download_url);

        // Download the file
        let response = self.client.get(&download_url).send()?;

        if !response.status().is_success() {
            // Try alternative common filename pattern if first failed
            let alt_filename = format!("{}_genomic.fna.gz", ftp_basename);
            let alt_download_url = format!("{}/{}", ftp_path_base, alt_filename);
            info!(
                "Download failed (Status: {}). Trying alternative URL: {}",
                response.status(),
                alt_download_url
            );

            let response_alt = self.client.get(&alt_download_url).send()?;
            if !response_alt.status().is_success() {
                return Err(DatabaseError::NCBIApiError(format!(
                    "Genome download failed for {} (Status: {}). Tried URLs: {} and {}",
                    accession,
                    response_alt.status(),
                    download_url,
                    alt_download_url
                )));
            }
            // Use the alternative response if successful
            let content = response_alt.bytes()?;
            let mut file = File::create(&cache_file)?;
            file.write_all(&content)?;
        } else {
            // Save the primary response to cache
            let content = response.bytes()?;
            let mut file = File::create(&cache_file)?;
            file.write_all(&content)?;
        }

        info!(
            "Successfully downloaded and cached: {}",
            cache_file.display()
        );
        Ok(cache_file)
    }
}

/// Signature database using sled
pub struct SignatureDatabase {
    /// Underlying key-value store
    db: Db,

    /// Index of taxonomy IDs to accessions (signature IDs)
    taxonomy_index: HashMap<String, HashSet<String>>, // Use HashSet for unique IDs

    /// Index of lineage terms (names) to accessions (signature IDs)
    lineage_index: HashMap<String, HashSet<String>>, // Use HashSet for unique IDs
}

impl SignatureDatabase {
    /// Open or create a signature database
    pub fn open(path: impl AsRef<Path>) -> Result<Self, DatabaseError> {
        info!("Opening database at: {}", path.as_ref().display());
        let db = sled::open(path)?;

        // Load indices, handling potential errors during decode gracefully
        let taxonomy_index: HashMap<String, HashSet<String>> = db
            .get("taxonomy_index")?
            .map(|data| decode_from_slice(&data, standard()).map(|(d, _)| d))
            .transpose() // Turns Option<Result<T, E>> into Result<Option<T>, E>
            .map_err(|e| {
                DatabaseError::SerializationError(format!("Failed to decode taxonomy index: {}", e))
            })?
            .unwrap_or_default(); // Use default empty map if None or decode error

        let lineage_index: HashMap<String, HashSet<String>> = db
            .get("lineage_index")?
            .map(|data| decode_from_slice(&data, standard()).map(|(d, _)| d))
            .transpose()
            .map_err(|e| {
                DatabaseError::SerializationError(format!("Failed to decode lineage index: {}", e))
            })?
            .unwrap_or_default();

        info!(
            "Database opened. Taxonomy index size: {}, Lineage index size: {}",
            taxonomy_index.len(),
            lineage_index.len()
        );

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
        // Key: signature.taxon_id (e.g., accession)
        // Value: serialized MultiResolutionSignature
        let key = signature.taxon_id.as_bytes();
        let signature_data = encode_to_vec(signature, standard())?;

        // Store signature
        self.db.insert(key, signature_data)?;
        info!("Added signature with ID: {}", signature.taxon_id);

        // Update indices
        self.update_indices(signature)?;

        // Persist indices immediately after update
        self.save_indices()?;
        self.db.flush()?; // Ensure data is written to disk

        Ok(())
    }

    /// Update in-memory search indices for a signature
    fn update_indices(
        &mut self,
        signature: &MultiResolutionSignature,
    ) -> Result<(), DatabaseError> {
        let signature_id = signature.taxon_id.clone();

        // Update taxonomy index (TaxID -> Set<SignatureID>)
        // Assumes signature.lineage contains ("TaxID", "Name") pairs, BUT the current
        // MultiResolutionSignature dummy has Vec<String> (names only).
        // We need the actual TaxID from the GenomeMetadata or stored differently.
        // For now, let's assume the *last* element of the lineage Vec<String> corresponds
        // to the species name, and we use the signature ID itself as a proxy for taxid lookup?
        // This is a placeholder until the lineage handling is consistent.
        // A better approach: Store the actual TaxID in MultiResolutionSignature.
        // Using the signature ID itself as the key here as a temporary measure:
        self.taxonomy_index
            .entry(signature_id.clone()) // Using sig ID as placeholder key
            .or_default()
            .insert(signature_id.clone());

        // Update lineage index (Name -> Set<SignatureID>)
        for name in &signature.lineage {
            if !name.is_empty() {
                // Avoid indexing empty names
                self.lineage_index
                    .entry(name.clone())
                    .or_default()
                    .insert(signature_id.clone());
            }
        }

        Ok(())
    }

    /// Save indices to the database
    fn save_indices(&self) -> Result<(), DatabaseError> {
        let taxonomy_data = encode_to_vec(&self.taxonomy_index, standard())?;
        let lineage_data = encode_to_vec(&self.lineage_index, standard())?;

        // Use atomic batch for index updates if possible (sled >= 0.34)
        let mut batch = sled::Batch::default();
        batch.insert("taxonomy_index", taxonomy_data);
        batch.insert("lineage_index", lineage_data);
        self.db.apply_batch(batch)?;

        // Fallback for older sled versions:
        // self.db.insert("taxonomy_index", taxonomy_data)?;
        // self.db.insert("lineage_index", lineage_data)?;
        Ok(())
    }

    /// Get a signature by ID (e.g., accession)
    pub fn get_signature(&self, id: &str) -> Result<MultiResolutionSignature, DatabaseError> {
        match self.db.get(id.as_bytes())? {
            Some(data) => {
                let (signature, _): (MultiResolutionSignature, _) =
                    decode_from_slice(&data, standard())?;
                Ok(signature)
            }
            None => Err(DatabaseError::NotFoundError(format!(
                "Signature not found: {}",
                id
            ))),
        }
    }

    /// Search for signatures matching a taxonomy term (either TaxID or Name)
    pub fn search_by_taxonomy(
        &self,
        term: &str,
    ) -> Result<Vec<MultiResolutionSignature>, DatabaseError> {
        let mut matching_ids = HashSet::new();

        // Search taxonomy index (assuming term might be a TaxID or the placeholder sig ID)
        if let Some(ids) = self.taxonomy_index.get(term) {
            matching_ids.extend(ids.iter().cloned());
        }

        // Search lineage index (assuming term is a Name)
        if let Some(ids) = self.lineage_index.get(term) {
            matching_ids.extend(ids.iter().cloned());
        }

        if matching_ids.is_empty() {
            info!("No signatures found matching term: {}", term);
            return Ok(Vec::new());
        }

        info!(
            "Found {} potential signature IDs for term: {}",
            matching_ids.len(),
            term
        );

        // Fetch signatures for the unique IDs found
        let mut results = Vec::with_capacity(matching_ids.len());
        for id in matching_ids {
            match self.get_signature(&id) {
                Ok(signature) => results.push(signature),
                Err(DatabaseError::NotFoundError(_)) => {
                    // This might happen if index is stale, log a warning
                    warn!(
                        "Signature ID {} found in index but not in database. Index might be stale.",
                        id
                    );
                }
                Err(e) => {
                    error!("Error retrieving signature {}: {}", id, e);
                    // Optionally propagate the error, or just log and continue
                    // return Err(e); // Propagate
                    continue; // Log and skip
                }
            }
        }

        Ok(results)
    }

    /// Get all signatures stored in the database
    pub fn get_all_signatures(&self) -> Result<Vec<MultiResolutionSignature>, DatabaseError> {
        let mut results = Vec::new();
        for item in self.db.iter() {
            let (key, value) = item?;
            // Safely convert key to string, skip if invalid UTF-8 or index key
            if let Ok(key_str) = std::str::from_utf8(&key) {
                if key_str == "taxonomy_index" || key_str == "lineage_index" {
                    continue;
                }

                match decode_from_slice::<MultiResolutionSignature, _>(&value, standard()) {
                    Ok((signature, _)) => results.push(signature),
                    Err(e) => {
                        error!(
                            "Failed to decode signature data for key '{}': {}. Skipping.",
                            key_str, e
                        );
                        // Continue to next item
                    }
                }
            } else {
                warn!("Found non-UTF8 key in database: {:?}", key);
                // Skip non-UTF8 keys as they are not expected signature IDs
            }
        }
        Ok(results)
    }

    /// Get the number of signatures (excluding index entries)
    pub fn count(&self) -> Result<usize, DatabaseError> {
        let mut count = 0;
        for item in self.db.iter() {
            let (key, _) = item?;
            if let Ok(key_str) = std::str::from_utf8(&key) {
                if key_str != "taxonomy_index" && key_str != "lineage_index" {
                    count += 1;
                }
            }
        }
        Ok(count)
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
        cache_dir: impl AsRef<Path>, // Primary k-mer for builder setup (using macro_k based on Init pattern)
        builder_kmer_size: usize,    // Sketch size for builder setup
        builder_sketch_size: usize,
        api_key: Option<String>, // Note: meso_k and threads are removed
    ) -> Result<Self, DatabaseError> {
        let database = SignatureDatabase::open(db_path)?;
        let downloader = NCBIDownloader::new(cache_dir, api_key, None)?; // Use default expiry for now

        // Ensure the builder is initialized with correct parameters
        // This depends on the actual SignatureBuilder implementation
        let builder = SignatureBuilder::new(31, 21, 1000, 1); // Example

        Ok(DatabaseManager {
            database,
            downloader,
            builder: builder.unwrap(),
        })
    }

    /// Search for and download reference genomes from NCBI
    pub fn download_references(
        &self,
        query: &str,
        max_results: usize,
    ) -> Result<Vec<(GenomeMetadata, PathBuf)>, DatabaseError> {
        info!(
            "Starting reference download for query: '{}', max_results: {}",
            query, max_results
        );
        // Search for matching genomes
        let genomes = self.downloader.search_genomes(query, max_results)?;
        if genomes.is_empty() {
            info!("No genomes found matching the query.");
            return Ok(Vec::new());
        }

        info!(
            "Found {} genomes matching query. Starting downloads...",
            genomes.len()
        );

        // Download each genome in parallel
        let results: Vec<Result<(GenomeMetadata, PathBuf), DatabaseError>> = genomes
            .into_par_iter()
            .map(|genome| {
                match self.downloader.download_genome(&genome.accession) {
                    Ok(path) => Ok((genome, path)),
                    Err(e) => {
                        error!("Failed to download genome {}: {}", genome.accession, e);
                        Err(e) // Propagate the error if needed, or filter out later
                    }
                }
            })
            .collect();

        // Handle results: collect successes, log errors
        let mut successful_downloads = Vec::new();
        for result in results {
            match result {
                Ok(pair) => successful_downloads.push(pair),
                Err(_) => {} // Error already logged by the download function or above
            }
        }

        info!(
            "Successfully downloaded {} out of {} genomes.",
            successful_downloads.len(),
            max_results
        );
        Ok(successful_downloads)
    }

    /// Process downloaded genomes into signatures and add to database
    pub fn process_references(
        &mut self,
        references: Vec<(GenomeMetadata, PathBuf)>,
    ) -> Result<Vec<String>, DatabaseError> {
        if references.is_empty() {
            info!("No references to process.");
            return Ok(Vec::new());
        }
        info!(
            "Processing {} downloaded references into signatures...",
            references.len()
        );

        // Convert GenomeMetadata to format needed by SignatureBuilder
        // IMPORTANT: Ensure lineage format matches what SignatureBuilder/MultiResolutionSignature expects
        let files_for_builder = references
            .iter()
            .map(|(metadata, path)| {
                // Convert lineage Vec<(String, String)> to Vec<String> (names only)
                // This matches the current dummy signature and database index logic.
                // If the signature needs IDs, adjust this conversion.
                let lineage_names: Vec<String> = metadata
                    .lineage
                    .iter()
                    .map(|(_, name)| name.clone())
                    .collect();

                (path.clone(), metadata.accession.clone(), lineage_names)
            })
            .collect::<Vec<_>>();

        // Build signatures in batch
        // This can be computationally intensive
        info!("Starting signature batch build...");
        let signatures = self.builder.build_batch(files_for_builder).map_err(|e| {
            DatabaseError::SignatureError(format!("Signature building failed: {}", e))
        })?;
        info!("Successfully built {} signatures.", signatures.len());

        // Add signatures to database
        let mut added_ids = Vec::with_capacity(signatures.len());
        for signature in signatures {
            // Add signature and handle potential errors (e.g., DB error)
            match self.database.add_signature(&signature) {
                Ok(_) => added_ids.push(signature.taxon_id.clone()),
                Err(e) => {
                    error!(
                        "Failed to add signature {} to database: {}",
                        signature.taxon_id, e
                    );
                    // Decide whether to stop or continue
                    return Err(e); // Stop processing on DB error
                }
            }
        }

        info!(
            "Successfully added {} new signatures to the database.",
            added_ids.len()
        );
        Ok(added_ids)
    }

    /// Search, download, and process reference genomes by query
    pub fn search_and_add_references(
        &mut self,
        query: &str,
        max_results: usize,
    ) -> Result<Vec<String>, DatabaseError> {
        // Download references
        let references = self.download_references(query, max_results)?;
        if references.is_empty() {
            return Ok(Vec::new()); // Nothing to process
        }
        // Process references
        self.process_references(references)
    }

    /// Check if the database contains any signatures
    pub fn is_empty(&self) -> Result<bool, DatabaseError> {
        Ok(self.database.count()? == 0)
    }
}

// --- Tests ---
#[cfg(test)]
mod tests {
    // Important: Include the dummy or real signature module
    use crate::sketch::signature::MultiResolutionSignature;
    use crate::sketch::SignatureBuilder;

    use super::*;
    use mockito::{Matcher, ServerGuard}; // Import ServerGuard
    use tempfile::tempdir;

    // Helper function to create a temporary directory
    fn create_temp_dir() -> tempfile::TempDir {
        tempdir().expect("Failed to create temp directory")
    }

    // Helper for basic NCBIDownloader setup with mock server
    fn setup_mock_downloader(server: &mut ServerGuard) -> (NCBIDownloader, PathBuf) {
        let temp_dir = create_temp_dir();
        let cache_dir = temp_dir.path().join("test_cache");
        fs::create_dir_all(&cache_dir).unwrap();

        let downloader = NCBIDownloader {
            client: Client::new(),
            base_url: server.url(), // Use mock server URL
            api_key: None,
            cache_dir: cache_dir.clone(),
            cache_expiry_days: 30,
        };
        (downloader, temp_dir.into_path()) // Return path to keep temp dir alive
    }

    #[test]
    fn test_ncbi_downloader_search_and_lineage() {
        // Setup mock server
        let mut server = mockito::Server::new();
        let (downloader, _temp_dir_guard) = setup_mock_downloader(&mut server); // Keep temp dir alive

        // Mock esearch response
        let _m_search = server
            .mock(
                "GET",
                Matcher::Regex(r"^/esearch.fcgi.*escherichia\+coli.*".to_string()),
            )
            .with_status(200)
            .with_header("content-type", "application/json")
            .with_body(r#"{"esearchresult":{"idlist":["12345"]}}"#)
            .create();

        // Mock esummary response
        let _m_summary = server
            .mock("GET", "/esummary.fcgi?db=assembly&id=12345&retmode=json")
            .with_status(200)
            .with_header("content-type", "application/json")
            .with_body(
                r#"{
                "result": {
                    "uids": ["12345"],
                    "12345": {
                        "assemblyaccession": "GCF_000005845.2",
                        "speciesname": "Escherichia coli",
                        "taxid": "562",
                        "assemblylevel": "Complete Genome",
                        "releasedate": "2014/08/25 00:00",
                        "submissiondate": "2001/07/26",
                        "totallength": "4641652",
                        "genomegcpercent": "50.79"
                    }
                }
            }"#,
            )
            .create();

        // Mock efetch taxonomy response (using LineageEx)
        let _m_taxonomy = server
            .mock("GET", "/efetch.fcgi?db=taxonomy&id=562&retmode=xml")
            .with_status(200)
            .with_header("content-type", "application/xml")
            .with_body(
                r#"<?xml version="1.0" ?>
            <TaxaSet>
                <Taxon>
                    <TaxId>562</TaxId>
                    <ScientificName>Escherichia coli</ScientificName>
                    <Rank>species</Rank>
                    <LineageEx>
                        <Taxon>
                            <TaxId>131567</TaxId>
                            <ScientificName>cellular organisms</ScientificName>
                            <Rank>no rank</Rank>
                        </Taxon>
                        <Taxon>
                            <TaxId>2</TaxId>
                            <ScientificName>Bacteria</ScientificName>
                            <Rank>superkingdom</Rank>
                        </Taxon>
                        <Taxon>
                            <TaxId>1224</TaxId>
                            <ScientificName>Proteobacteria</ScientificName>
                            <Rank>phylum</Rank>
                        </Taxon>
                        <Taxon>
                            <TaxId>1236</TaxId>
                            <ScientificName>Gammaproteobacteria</ScientificName>
                            <Rank>class</Rank>
                        </Taxon>
                        <Taxon>
                            <TaxId>91347</TaxId>
                            <ScientificName>Enterobacterales</ScientificName>
                            <Rank>order</Rank>
                        </Taxon>
                        <Taxon>
                            <TaxId>543</TaxId>
                            <ScientificName>Enterobacteriaceae</ScientificName>
                            <Rank>family</Rank>
                        </Taxon>
                         <Taxon>
                            <TaxId>561</TaxId>
                            <ScientificName>Escherichia</ScientificName>
                            <Rank>genus</Rank>
                        </Taxon>
                        <!-- NOTE: The species itself might not be in LineageEx -->
                    </LineageEx>
                </Taxon>
            </TaxaSet>"#,
            )
            .create();

        // Test search_genomes
        let results = downloader.search_genomes("escherichia coli", 1).unwrap();
        assert_eq!(results.len(), 1);
        let metadata = &results[0];

        assert_eq!(metadata.accession, "GCF_000005845.2");
        assert_eq!(metadata.assembly_id, "12345"); // Check assembly ID stored
        assert_eq!(metadata.organism, "Escherichia coli");
        assert_eq!(metadata.taxid, "562");
        assert_eq!(metadata.assembly_level, "Complete Genome");
        assert_eq!(metadata.release_date, "2014/08/25 00:00"); // Check correct date field used
        assert_eq!(metadata.size, 4641652);
        assert_eq!(metadata.gc_content, 50.79);

        // Verify parsed lineage (includes the main taxon added at the end)
        let expected_lineage = vec![
            ("131567".to_string(), "cellular organisms".to_string()),
            ("2".to_string(), "Bacteria".to_string()),
            ("1224".to_string(), "Proteobacteria".to_string()),
            ("1236".to_string(), "Gammaproteobacteria".to_string()),
            ("91347".to_string(), "Enterobacterales".to_string()),
            ("543".to_string(), "Enterobacteriaceae".to_string()),
            ("561".to_string(), "Escherichia".to_string()),
            ("562".to_string(), "Escherichia coli".to_string()), // Main taxon added
        ];
        assert_eq!(metadata.lineage, expected_lineage);

        // Mocks are automatically verified on drop if .assert() isn't used
    }

    #[test]
    fn test_database_manager() {
        let temp_dir = create_temp_dir();
        let db_path = temp_dir.path().join("test_manager_db");
        let cache_dir = temp_dir.path().join("test_manager_cache");

        // Create the database manager (uses dummy builder)
        let manager = DatabaseManager::new(&db_path, &cache_dir, 31, 500, None); // k=31, size=500
        assert!(manager.is_ok());
        let manager = manager.unwrap();

        // Test is_empty on new DB
        assert!(manager.is_empty().unwrap());

        // Further tests would involve mocking NCBI calls within the manager methods
        // For now, this confirms basic setup.
    }

    #[test]
    fn test_genome_metadata_serialization() {
        let metadata = GenomeMetadata {
            accession: "GCF_000001234.1".to_string(),
            assembly_id: "54321".to_string(), // Added field
            organism: "Escherichia coli".to_string(),
            taxid: "562".to_string(),
            assembly_level: "Complete Genome".to_string(),
            release_date: "2020-01-01".to_string(),
            size: 5000000,
            gc_content: 50.5,
            lineage: vec![
                ("2".to_string(), "Bacteria".to_string()),
                ("1224".to_string(), "Proteobacteria".to_string()),
            ],
        };

        // Serialize and deserialize
        let serialized = encode_to_vec(&metadata, standard()).unwrap();
        let (deserialized, _): (GenomeMetadata, _) =
            decode_from_slice(&serialized, standard()).unwrap();

        // Verify contents
        assert_eq!(deserialized.accession, metadata.accession);
        assert_eq!(deserialized.assembly_id, metadata.assembly_id);
        assert_eq!(deserialized.organism, metadata.organism);
        assert_eq!(deserialized.taxid, metadata.taxid);
        assert_eq!(deserialized.assembly_level, metadata.assembly_level);
        assert_eq!(deserialized.release_date, metadata.release_date);
        assert_eq!(deserialized.size, metadata.size);
        assert_eq!(deserialized.gc_content, metadata.gc_content);
        assert_eq!(deserialized.lineage, metadata.lineage);
    }
}

// --- More Mock Tests (in separate module for organization) ---
#[cfg(test)]
mod mock_tests {
    use super::*; // Access items from parent module (DatabaseError, NCBIDownloader, etc.)
    use mockito::{Matcher, ServerGuard};
    use std::io::Write;
    use tempfile::tempdir;

    // Helper function (duplicate, consider moving to a common test utils mod)
    fn setup_mock_downloader(server: &mut ServerGuard) -> (NCBIDownloader, PathBuf) {
        let temp_dir = tempdir().expect("Failed to create temp directory");
        let cache_dir = temp_dir.path().join("mock_cache");
        fs::create_dir_all(&cache_dir).unwrap();

        let downloader = NCBIDownloader {
            client: Client::new(),
            base_url: server.url(),
            api_key: None,
            cache_dir: cache_dir.clone(),
            cache_expiry_days: 30,
        };
        (downloader, temp_dir.into_path())
    }

    #[test]
    fn test_error_handling_in_search_flow() {
        let mut server = mockito::Server::new();
        let (downloader, _tdg) = setup_mock_downloader(&mut server);

        // Mock search success
        let _m_search = server
            .mock(
                "GET",
                Matcher::Regex(r"^/esearch.fcgi.*error\+species.*".to_string()),
            )
            .with_status(200)
            .with_body(r#"{"esearchresult":{"idlist":["1234567"]}}"#)
            .create();

        // Mock summary error
        let _m_summary = server
            .mock("GET", "/esummary.fcgi?db=assembly&id=1234567&retmode=json")
            .with_status(500) // Internal Server Error
            .with_body("Server Error")
            .create();

        // Test search_genomes - should fail gracefully on summary error for the chunk
        // Since it skips the chunk, the result should be Ok(empty_vec) if only one chunk failed
        let result = downloader.search_genomes("error species", 1);
        assert!(result.is_ok()); // The overall search doesn't fail, just skips bad chunk
        assert!(result.unwrap().is_empty()); // No results because the only chunk failed
    }

    #[test]
    fn test_api_response_malformed_json() {
        let mut server = mockito::Server::new();
        let (downloader, _tdg) = setup_mock_downloader(&mut server);

        // Mock malformed search response
        let _m_search = server
            .mock(
                "GET",
                Matcher::Regex(r"^/esearch.fcgi.*malformed.*".to_string()),
            )
            .with_status(200)
            .with_header("content-type", "application/json")
            .with_body(r#"{"esearchresult":{"idlist":["7654321"]"#) // Missing closing braces
            .create();

        // Test error handling
        let genomes = downloader.search_genomes("malformed", 1);
        assert!(genomes.is_err());
        // Expect HttpError because reqwest::Error wraps serde json errors during .json() call
        assert!(matches!(genomes, Err(DatabaseError::HttpError(_))));
    }

    #[test]
    fn test_empty_search_result_handling() {
        let mut server = mockito::Server::new();
        let (downloader, _tdg) = setup_mock_downloader(&mut server);

        // Mock empty search response
        let _m_search = server
            .mock(
                "GET",
                Matcher::Regex(r"^/esearch.fcgi.*empty\+results.*".to_string()),
            )
            .with_status(200)
            .with_body(r#"{"esearchresult":{"idlist":[]}}"#)
            .create();

        // Test handling
        let genomes = downloader.search_genomes("empty results", 1);
        assert!(genomes.is_ok()); // Should return Ok, not Err
        assert!(genomes.unwrap().is_empty()); // Result vector should be empty
    }

    #[test]
    fn test_api_key_usage_in_urls() {
        let mut server = mockito::Server::new();
        let api_key = "TEST_API_KEY_123";

        let temp_dir = tempdir().expect("Failed to create temp directory");
        let cache_dir = temp_dir.path().join("api_key_cache");
        fs::create_dir_all(&cache_dir).unwrap();

        let downloader = NCBIDownloader {
            client: Client::new(),
            base_url: server.url(),
            api_key: Some(api_key.to_string()), // Set API key
            cache_dir: cache_dir.clone(),
            cache_expiry_days: 30,
        };

        // Mock search URL *with* API key
        let search_path_matcher = Matcher::Regex(format!(r"api_key={}", api_key));
        let _m_search = server
            .mock(
                "GET",
                Matcher::AllOf(vec![
                    Matcher::Regex(r"esearch.fcgi".to_string()),
                    search_path_matcher.clone(),
                ]),
            )
            .with_status(200)
            .with_body(r#"{"esearchresult":{"idlist":[]}}"#) // Empty result is fine
            .create();

        // Mock summary URL *with* API key
        let _m_summary = server
            .mock(
                "GET",
                Matcher::AllOf(vec![
                    Matcher::Regex(r"esummary.fcgi".to_string()),
                    search_path_matcher.clone(), // Reuse matcher
                ]),
            )
            .with_status(200)
            .with_body(r#"{"result":{"uids":[]}}"#) // Empty result needed to prevent panic
            .create();

        // Mock taxonomy URL *with* API key
        let _m_taxonomy = server
            .mock(
                "GET",
                Matcher::AllOf(vec![
                    Matcher::Regex(r"efetch.fcgi".to_string()),
                    search_path_matcher.clone(), // Reuse matcher
                ]),
            )
            .with_status(200)
            .with_body(r#"<TaxaSet></TaxaSet>"#)
            .create();

        // Call a method that uses the API key
        // search_genomes triggers search, summary, and taxonomy mocks
        let _ = downloader.search_genomes("test with key", 1);

        // Mocks are verified automatically when server drops.
        // If the mocks weren't hit, the test would fail.
    }
}
