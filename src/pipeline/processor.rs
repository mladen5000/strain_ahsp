use crate::adaptive::classifier::{AdaptiveClassifier, Classification, TaxonomicLevel};
use crate::database::DatabaseManager;
// Correct signature imports (assuming these are the canonical types)
use crate::sketch::signature::{KmerSignature, MultiResolutionSignature}; // Add VariantProfile if needed by MultiResSig
use log::{debug, error, info, warn}; // Include debug
use needletail::parse_fastx_file; // Correct import
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::{self, File}; // Use fs for create_dir_all
use std::io::{self, BufWriter};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::time::Instant;
use thiserror::Error;

// --- Error Type ---
#[derive(Error, Debug)]
pub enum ProcessingError {
    #[error("IO error: {0}")]
    IoError(#[from] io::Error),

    #[error("FASTQ parsing error: {0}")]
    FastqError(String),

    // Assuming SignatureError might come from sketch module or local operations
    #[error("Signature error: {0}")]
    SignatureError(String),

    #[error("Classification error: {0}")]
    ClassificationError(String),

    #[error("Strain estimation error: {0}")]
    StrainEstimationError(String),

    #[error("Database error: {0}")]
    DatabaseError(String),

    #[error("Needletail parsing error: {0}")]
    NeedletailError(#[from] needletail::errors::ParseError),

    #[error("Classifier not initialized")] // Added specific error
    ClassifierUninitialized,

    #[error("Failed to lock mutex: {0}")] // Added specific error
    MutexLockError(String),

    #[error("No valid reference signatures could be loaded or converted")] // Added specific error
    NoValidReferences,
}

// Custom conversion from Mutex PoisonError
impl<T> From<std::sync::PoisonError<T>> for ProcessingError {
    fn from(err: std::sync::PoisonError<T>) -> Self {
        ProcessingError::MutexLockError(err.to_string())
    }
}

// --- Structs (QC Params, Metrics, Results) ---
// (Keeping definitions as provided)

/// Quality control parameters for FASTQ processing
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityControlParams {
    pub min_avg_quality: f64,
    pub min_length: usize,
    pub trim_quality: u8,
    pub max_n_percent: f64,
}

impl Default for QualityControlParams {
    fn default() -> Self {
        QualityControlParams {
            min_avg_quality: 20.0,
            min_length: 50,
            trim_quality: 15,
            max_n_percent: 5.0,
        }
    }
}

/// Processing metrics
#[derive(Debug, Clone, Serialize, Deserialize, Default)] // Added Default
pub struct ProcessingMetrics {
    pub total_reads: usize,
    pub passed_reads: usize,
    pub total_bases: usize,
    pub passed_bases: usize,
    pub avg_read_length: f64,
    pub processing_time_seconds: f64,
}

/// Sample classification results
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ClassificationResults {
    pub sample_id: String,
    pub metrics: ProcessingMetrics,
    pub classifications: Vec<Classification>, // Store top classification(s)
    pub strain_abundances: HashMap<String, (f64, f64)>,
    pub results_file: Option<PathBuf>,
}

// --- FastqProcessor ---

/// FASTQ processing pipeline
pub struct FastqProcessor {
    qc_params: QualityControlParams,
    threads: usize,
    chunk_size: usize,
    macro_k: usize,
    meso_k: usize,
    sketch_size: usize,
    db_manager: DatabaseManager,
    // Store the classifier directly, potentially wrapped in Arc if needed elsewhere
    classifier: Option<AdaptiveClassifier>,
}

impl FastqProcessor {
    /// Create a new FASTQ processor instance.
    pub fn new(
        db_path: impl AsRef<Path>,
        cache_dir: impl AsRef<Path>,
        threads: usize,
        macro_k: usize,
        meso_k: usize,
        sketch_size: usize,
        qc_params: Option<QualityControlParams>,
        api_key: Option<String>,
    ) -> Result<Self, ProcessingError> {
        info!("Initializing FastqProcessor...");
        // Basic validation
        if macro_k < 3 || meso_k < 3 || sketch_size == 0 {
            return Err(ProcessingError::SignatureError(format!(
                 "Invalid parameters: macro_k={}, meso_k={}, sketch_size={}. K must be >= 3, sketch > 0.",
                 macro_k, meso_k, sketch_size
            )));
        }

        let db_manager = DatabaseManager::new(
            db_path.as_ref(), // Use as_ref()
            cache_dir.as_ref(),
            sketch_size,
            threads,
            api_key,
        )
        .map_err(|e| ProcessingError::DatabaseError(format!("DB Manager init failed: {}", e)))?;
        info!("DatabaseManager initialized.");

        Ok(FastqProcessor {
            qc_params: qc_params.unwrap_or_default(),
            threads: threads.max(1), // Ensure at least one thread
            chunk_size: 100_000,     // Default chunk size
            macro_k,
            meso_k,
            sketch_size,
            db_manager,
            classifier: None, // Initialize classifier separately
        })
    }

    /// Initialize the classifier by loading and preparing reference signatures.
    /// Must be called after `new` and before `process_file`.
    pub fn init_classifier(&mut self) -> Result<(), ProcessingError> {
        if self.classifier.is_some() {
            info!("Classifier already initialized.");
            return Ok(());
        }
        info!("Initializing classifier: Fetching reference signatures...");

        // Assume get_all_signatures returns Vec<Arc<SomeDbType>>
        // **CRITICAL:** The structure of `db_sig_arc` MUST match expectations below.
        let db_references = self.db_manager.database.get_all_signatures().map_err(|e| {
            ProcessingError::DatabaseError(format!("Failed to get signatures: {}", e))
        })?;

        if db_references.is_empty() {
            return Err(ProcessingError::NoValidReferences);
        }
        info!(
            "Loaded {} reference signatures from DB.",
            db_references.len()
        );

        // Convert DB signatures to the sketch::signature::MultiResolutionSignature type
        let mut sketch_signatures = Vec::with_capacity(db_references.len());
        let mut skipped_count = 0;
        for db_sig_arc in db_references {
            // --- Conversion Logic ---
            // This is the most fragile part. Assumes `db_sig_arc` has the necessary fields.
            // Ensure `KmerSignature` and `Option<VariantProfile>` types match.
            if let (Some(macro_sig), Some(meso_sig)) = (
                db_sig_arc.macro_signature.as_ref(), // Borrow instead of cloning immediately
                db_sig_arc.meso_signature.as_ref(),
            ) {
                // Construct the required signature type
                let sketch_sig = MultiResolutionSignature {
                    taxon_id: db_sig_arc.taxon_id.clone(),
                    lineage: db_sig_arc.lineage.clone(),
                    macro_signature: macro_sig.clone(), // Now clone
                    meso_signature: meso_sig.clone(),   // Now clone
                    // Clone micro if present, otherwise None
                    micro_signature: db_sig_arc.micro_signature.clone(),
                    // Initialize weights or other fields as needed by the struct definition
                    // Assuming 'weights' is not a field based on previous errors,
                    // If a field like 'levels' exists, initialize it appropriately e.g., levels: vec![]
                    weights: HashMap::new(), // Or remove if field doesn't exist
                    levels: vec![],          // Use if 'levels' is the field
                };
                sketch_signatures.push(Arc::new(sketch_sig));
            } else {
                warn!(
                    "Skipping reference signature '{}' due to missing macro or meso signature.",
                    db_sig_arc.taxon_id
                );
                skipped_count += 1;
            }
            // --- End Conversion Logic ---
        }

        if sketch_signatures.is_empty() {
            error!(
                "Failed to convert any DB signatures. {} signatures were skipped.",
                skipped_count
            );
            return Err(ProcessingError::NoValidReferences);
        }
        info!(
            "Successfully converted {} signatures ({} skipped). Creating classifier...",
            sketch_signatures.len(),
            skipped_count
        );

        // Create the classifier instance
        let classifier = AdaptiveClassifier::new(sketch_signatures, None, None) // Pass Vec<Arc<MultiResolutionSignature>>
            .map_err(|e| {
                ProcessingError::ClassificationError(format!("Classifier creation failed: {}", e))
            })?;
        info!("Classifier initialized successfully.");

        self.classifier = Some(classifier);
        Ok(())
    }

    /// Get a reference to the initialized classifier.
    /// Returns an error if the classifier has not been initialized.
    fn get_classifier(&self) -> Result<&AdaptiveClassifier, ProcessingError> {
        self.classifier
            .as_ref()
            .ok_or(ProcessingError::ClassifierUninitialized)
    }

    /// Process a single FASTQ file.
    pub fn process_file(
        &self,
        fastq_path: impl AsRef<Path>,
        sample_id: &str,
        output_dir: impl AsRef<Path>,
    ) -> Result<ClassificationResults, ProcessingError> {
        let start_time = Instant::now();
        info!(
            "Starting processing for sample '{}' from file: {}",
            sample_id,
            fastq_path.as_ref().display()
        );

        let classifier = self.get_classifier()?; // Ensure classifier is ready

        // Prepare output directory
        let output_path = output_dir.as_ref();
        fs::create_dir_all(output_path)?;
        info!("Output directory set to: {}", output_path.display());

        // Initialize metrics and signature in Arcs for thread safety
        let metrics = Arc::new(Mutex::new(ProcessingMetrics::default()));
        let sample_signature = self.create_empty_sample_signature(sample_id)?;

        // Read and process the file in chunks
        self.read_and_process_reads(fastq_path.as_ref(), &metrics, &sample_signature)?;

        let elapsed_read_time = start_time.elapsed();
        info!(
            "Read processing completed in {:.2} seconds.",
            elapsed_read_time.as_secs_f64()
        );

        // Finalize metrics
        let final_metrics = self.finalize_metrics(metrics, elapsed_read_time.as_secs_f64())?;

        // Get the final signature
        let final_signature = Arc::try_unwrap(sample_signature) // Try to get ownership
            .map_err(|_| ProcessingError::MutexLockError("Failed to unwrap Arc".to_string()))?
            .into_inner()?; // Get ownership from Mutex

        // Classify the signature
        info!("Classifying final sample signature...");
        // Assuming classify returns the single best result
        let classification = classifier.classify(&final_signature)?;
        let classifications = vec![classification]; // Wrap in Vec for consistency

        // Estimate strain abundances
        let strain_abundances = self.run_strain_estimation(&final_signature, &classifications)?;

        // Assemble and save results
        let results = self.assemble_and_save_results(
            sample_id,
            final_metrics,
            classifications,
            strain_abundances,
            output_path,
        )?;

        let total_time = start_time.elapsed().as_secs_f64();
        info!(
            "Total processing for sample '{}' finished in {:.2} seconds.",
            sample_id, total_time
        );
        self.log_summary_metrics(&results.metrics);

        Ok(results)
    }

    /// Creates an empty MultiResolutionSignature for a new sample.
    fn create_empty_sample_signature(
        &self,
        sample_id: &str,
    ) -> Result<Arc<Mutex<MultiResolutionSignature>>, ProcessingError> {
        // Initialize KmerSignatures using struct literals if ::new is unavailable
        let macro_sig = KmerSignature {
            kmers: Vec::with_capacity(self.sketch_size),
            counts: Vec::with_capacity(self.sketch_size),
            total_kmers: 0,
        };
        let meso_sig = KmerSignature {
            kmers: Vec::with_capacity(self.sketch_size),
            counts: Vec::with_capacity(self.sketch_size),
            total_kmers: 0,
        };

        // Create the MultiResolutionSignature using struct literal
        let signature = MultiResolutionSignature {
            taxon_id: sample_id.to_string(),
            lineage: Vec::new(), // Sample lineage is unknown initially
            macro_signature: macro_sig,
            meso_signature: meso_sig,
            micro_signature: None, // Initialize as None
            weights: HashMap::new(), // Or remove if field doesn't exist
                                   // levels: vec![], // Use if 'levels' is the field
        };

        Ok(Arc::new(Mutex::new(signature)))
    }

    /// Reads FASTQ records and processes them in parallel chunks.
    fn read_and_process_reads(
        &self,
        fastq_path: &Path,
        metrics: &Arc<Mutex<ProcessingMetrics>>,
        signature: &Arc<Mutex<MultiResolutionSignature>>,
    ) -> Result<(), ProcessingError> {
        let mut reader = parse_fastx_file(fastq_path)?;
        let mut current_chunk = Vec::with_capacity(self.chunk_size);

        while let Some(record_result) = reader.next() {
            let record = record_result?;
            current_chunk.push((record.seq().to_vec(), record.qual().map(|q| q.to_vec())));

            if current_chunk.len() >= self.chunk_size {
                self.process_chunk(t_chunk, metrics, signature)?;
                current_chunk.clear();
            }
        }

        if !current_chunk.is_empty() {
            self.process_chunk(t_chunk, metrics, signature)?;
        }
        Ok(())
    }

    /// Process a chunk of reads in parallel: apply QC and update the shared signature.
    fn process_chunk(
        &self,
        chunk: &[(Vec<u8>, Option<Vec<u8>>)],
        metrics: &Arc<Mutex<ProcessingMetrics>>,
        signature: &Arc<Mutex<MultiResolutionSignature>>,
    ) -> Result<(), ProcessingError> {
        // Use Rayon's parallel bridge for iterators combined with try_for_each for error handling
        chunk
            .par_iter()
            .try_for_each(|(seq, qual)| -> Result<(), ProcessingError> {
                // Ensure closure returns Result
                let original_len = seq.len();
                let processed_seq_opt = self.apply_quality_control(seq, qual.as_ref());

                // --- Metrics Update ---
                let passed_len = processed_seq_opt.as_ref().map_or(0, |v| v.len());
                let passed_qc = processed_seq_opt.is_some();
                {
                    // Scoped lock for metrics
                    let mut metrics_guard = metrics.lock()?; // Use '?' for PoisonError conversion
                    metrics_guard.total_reads += 1;
                    metrics_guard.total_bases += original_len;
                    if passed_qc {
                        metrics_guard.passed_reads += 1;
                        metrics_guard.passed_bases += passed_len;
                    }
                } // Mutex guard dropped

                // --- Signature Update ---
                if let Some(processed_seq) = processed_seq_opt {
                    if !processed_seq.is_empty() {
                        let mut sig_guard = signature.lock()?; // Use '?' for PoisonError conversion
                                                               // Call add_sequence on the inner KmerSignature fields
                                                               // Assuming KmerSignature has add_sequence method returning Result<(), SigError>
                        sig_guard
                            .macro_signature
                            .add_sequence(&processed_seq)
                            .map_err(|e| {
                                ProcessingError::SignatureError(format!(
                                    "Macro sig update failed: {}",
                                    e
                                ))
                            })?;
                        sig_guard
                            .meso_signature
                            .add_sequence(&processed_seq)
                            .map_err(|e| {
                                ProcessingError::SignatureError(format!(
                                    "Meso sig update failed: {}",
                                    e
                                ))
                            })?;
                    }
                }
                Ok(()) // Success for this item
            }) // try_for_each propagates first Err
    }

    /// Apply quality control filters to a single read.
    fn apply_quality_control(&self, seq: &[u8], qual: Option<&Vec<u8>>) -> Option<Vec<u8>> {
        // Apply QC logic (keeping previous implementation as it seems reasonable)
        if seq.len() < self.qc_params.min_length {
            return None;
        }

        let n_count = seq
            .iter()
            .filter(|&&base| base == b'N' || base == b'n')
            .count();
        let n_percent = 100.0 * n_count as f64 / seq.len() as f64;
        if n_percent > self.qc_params.max_n_percent {
            return None;
        }

        if let Some(qual_vec) = qual {
            if qual_vec.len() != seq.len() {
                error!(
                    "Seq/Qual length mismatch ({} vs {}). Discarding.",
                    seq.len(),
                    qual_vec.len()
                );
                return None;
            }
            if qual_vec.is_empty() {
                return None;
            }

            let avg_quality = qual_vec
                .iter()
                .map(|&q| q.saturating_sub(33) as f64)
                .sum::<f64>()
                / qual_vec.len() as f64;
            if avg_quality < self.qc_params.min_avg_quality {
                return None;
            }

            let mut trim_start = 0;
            let mut trim_end = seq.len();
            let mut found_start = false;
            for (i, &q) in qual_vec.iter().enumerate() {
                if q.saturating_sub(33) >= self.qc_params.trim_quality {
                    trim_start = i;
                    found_start = true;
                    break;
                }
            }
            if !found_start {
                return None;
            }

            let mut found_end = false;
            for i in (trim_start..qual_vec.len()).rev() {
                if qual_vec[i].saturating_sub(33) >= self.qc_params.trim_quality {
                    trim_end = i + 1;
                    found_end = true;
                    break;
                }
            }
            if !found_end {
                return None;
            }

            if trim_start >= trim_end || (trim_end - trim_start) < self.qc_params.min_length {
                return None;
            }

            Some(seq[trim_start..trim_end].to_vec())
        } else {
            Some(seq.to_vec()) // Passed length/N%, no quality scores
        }
    }

    /// Finalize metrics calculations (processing time, average length).
    fn finalize_metrics(
        &self,
        metrics_arc: Arc<Mutex<ProcessingMetrics>>,
        processing_time_seconds: f64,
    ) -> Result<ProcessingMetrics, ProcessingError> {
        let mut metrics_guard = metrics_arc.lock()?; // Use '?'
        metrics_guard.processing_time_seconds = processing_time_seconds;
        if metrics_guard.passed_reads > 0 {
            metrics_guard.avg_read_length =
                metrics_guard.passed_bases as f64 / metrics_guard.passed_reads as f64;
        } else {
            metrics_guard.avg_read_length = 0.0; // Avoid division by zero
        }
        // Clone the final metrics struct
        Ok(metrics_guard.clone())
    }

    /// Get hierarchical classifications (currently just the best hit).
    fn get_hierarchical_classifications(
        &self,
        signature: &MultiResolutionSignature,
        classifier: &AdaptiveClassifier,
    ) -> Result<Vec<Classification>, ProcessingError> {
        let best_classification = classifier.classify(signature)?; // Propagate error
        Ok(vec![best_classification])
    }

    /// Estimate strain abundances if classification is specific enough.
    fn run_strain_estimation(
        &self,
        final_signature: &MultiResolutionSignature,
        classifications: &[Classification],
    ) -> Result<HashMap<String, (f64, f64)>, ProcessingError> {
        if let Some(cls) = classifications.first() {
            if cls.level <= TaxonomicLevel::Species {
                info!(
                    "Classification at/below species level ('{}', {:?}). Estimating strains...",
                    cls.taxon_id, cls.level
                );
                let classifier = self.get_classifier()?; // Get classifier ref again
                self.estimate_strain_abundances(final_signature, classifier, &cls.taxon_id)
            } else {
                info!(
                    "Classification level ({:?}) is above Species. Skipping strain estimation.",
                    cls.level
                );
                Ok(HashMap::new())
            }
        } else {
            warn!("No classification available to guide strain estimation.");
            Ok(HashMap::new())
        }
    }

    /// Estimate relative abundances of strains related to the classified species.
    fn estimate_strain_abundances(
        &self,
        sample_signature: &MultiResolutionSignature,
        classifier: &AdaptiveClassifier,
        target_species_id: &str,
    ) -> Result<HashMap<String, (f64, f64)>, ProcessingError> {
        debug!(
            "Filtering references for strain estimation relative to {}",
            target_species_id
        );

        // Filter references based on lineage and ensure they are not the target species itself
        let relevant_strains: Vec<_> = classifier
            .references // Assuming Vec<Arc<MultiResolutionSignature>>
            .iter()
            .filter(|ref_sig| {
                ref_sig.lineage.contains(&target_species_id.to_string())
                    && ref_sig.taxon_id != target_species_id
            })
            .collect();

        if relevant_strains.is_empty() {
            warn!(
                "No potential reference strains found downstream of target {}.",
                target_species_id
            );
            return Ok(HashMap::new());
        }
        info!(
            "Found {} potential reference strains for target {}.",
            relevant_strains.len(),
            target_species_id
        );

        // Calculate similarities
        let mut similarities = HashMap::new();
        let mut total_similarity = 0.0;
        for strain_sig in relevant_strains {
            // Use overall weighted similarity (None uses weights defined in signature)
            // Or specify a level like Some(ResolutionLevel::Meso) if desired
            let sim = sample_signature.similarity(strain_sig, None);
            if sim > f64::EPSILON {
                // Use epsilon for float comparison
                // Clone the taxon_id string for the HashMap key
                similarities.insert(strain_sig.taxon_id.clone(), sim);
                total_similarity += sim;
            }
        }

        // Normalize to get abundances
        let mut abundances = HashMap::new();
        if total_similarity > f64::EPSILON {
            for (id, sim) in similarities {
                // `id` is already String here
                let abundance = sim / total_similarity;
                let confidence = sim.min(1.0); // Use similarity as a proxy for confidence (0-1) - very basic!
                abundances.insert(id, (abundance, confidence)); // id is moved into map here
            }
        } else {
            info!("Total similarity to relevant strains is negligible.");
        }

        // Sort abundances for logging (optional)
        let mut sorted_abundances: Vec<_> = abundances.iter().collect();
        sorted_abundances.sort_by(|a, b| {
            b.1 .0
                .partial_cmp(&a.1 .0)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        for (id, (abundance, confidence)) in sorted_abundances {
            info!(
                "  Strain {}: Relative Abundance ~{:.2}%, Confidence {:.2}",
                id,
                abundance * 100.0,
                confidence // Confidence is now 0-1
            );
        }

        Ok(abundances) // Return the calculated abundances map
    }

    /// Assemble the final results structure and save it to a JSON file.
    fn assemble_and_save_results(
        &self,
        sample_id: &str,
        metrics: ProcessingMetrics,
        classifications: Vec<Classification>,
        strain_abundances: HashMap<String, (f64, f64)>,
        output_path: &Path,
    ) -> Result<ClassificationResults, ProcessingError> {
        let results_file_path = output_path.join(format!("{}_results.json", sample_id));

        let results = ClassificationResults {
            sample_id: sample_id.to_string(),
            metrics,
            classifications,
            strain_abundances,
            results_file: Some(results_file_path.clone()),
        };

        info!("Writing final results to {}", results_file_path.display());
        let file = File::create(&results_file_path)?;
        let writer = BufWriter::new(file);
        // Use serde_json to write pretty-printed JSON
        serde_json::to_writer_pretty(writer, &results)
            .map_err(|e| ProcessingError::IoError(io::Error::new(io::ErrorKind::Other, e)))?;

        Ok(results)
    }

    /// Log summary metrics after processing.
    fn log_summary_metrics(&self, metrics: &ProcessingMetrics) {
        info!(
            "Summary Metrics: Reads Passed QC: {} / {} ({:.1}%)",
            metrics.passed_reads,
            metrics.total_reads,
            100.0 * metrics.passed_reads as f64 / metrics.total_reads.max(1) as f64
        );
        info!(
            "Summary Metrics: Avg Read Length (Passed QC): {:.1} bp",
            metrics.avg_read_length
        );
        info!(
            "Summary Metrics: Total Processing Time: {:.2} seconds",
            metrics.processing_time_seconds
        );
    }
} // end impl FastqProcessor

// --- Report Generation --- (Keeping definition as provided)

/// Generate a formatted text report from the classification results.
pub fn generate_report(results: &ClassificationResults) -> Result<String, ProcessingError> {
    let mut report = String::new();

    // Header
    report.push_str(&format!(
        "Classification Report for Sample: {}\n",
        results.sample_id
    ));
    report.push_str("=================================================\n\n");

    // Metrics Section
    report.push_str("Processing Metrics:\n");
    report.push_str(&format!(
        "  Total reads processed: {}\n",
        results.metrics.total_reads
    ));
    report.push_str(&format!(
        "  Reads passed QC: {} ({:.1}%)\n",
        results.metrics.passed_reads,
        100.0 * results.metrics.passed_reads as f64 / results.metrics.total_reads.max(1) as f64
    ));
    report.push_str(&format!(
        "  Bases passed QC: {}\n",
        results.metrics.passed_bases
    ));
    report.push_str(&format!(
        "  Average read length (passed QC): {:.1} bp\n",
        results.metrics.avg_read_length
    ));
    report.push_str(&format!(
        "  Processing time: {:.2} seconds\n\n",
        results.metrics.processing_time_seconds
    ));

    // Classification Section
    if results.classifications.is_empty() {
        report.push_str("Classification Results: No confident classification found.\n\n");
    } else {
        report.push_str("Classification Results (Top Hit):\n");
        if let Some(classification) = results.classifications.first() {
            report.push_str(&format!("  Taxon ID: {}\n", classification.taxon_id));
            report.push_str(&format!("  Taxonomic level: {:?}\n", classification.level));
            report.push_str(&format!("  Confidence: {:.4}\n", classification.confidence));

            if !classification.lineage.is_empty() {
                report.push_str("  Lineage: ");
                report.push_str(&classification.lineage.join(" > "));
                report.push('\n');
            } else {
                report.push_str("  Lineage: N/A\n");
            }

            report.push_str("  Similarity scores:\n");
            if classification.similarity_scores.is_empty() {
                report.push_str("    N/A\n");
            } else {
                // Sort scores by value (descending)
                let mut sorted_scores: Vec<_> = classification.similarity_scores.iter().collect();
                sorted_scores.sort_by(|(_, score_a), (_, score_b)| {
                    score_b
                        .partial_cmp(score_a)
                        .unwrap_or(std::cmp::Ordering::Equal)
                });
                for (level, score) in sorted_scores {
                    report.push_str(&format!("    {:?}: {:.4}\n", level, score));
                    // Assumes level is Debug
                }
            }
            report.push('\n');
        }
    }

    // Strain Abundance Section
    if !results.strain_abundances.is_empty() {
        report.push_str("Strain Abundance Estimates (relative within classified group):\n");
        let mut strains: Vec<_> = results.strain_abundances.iter().collect();
        strains.sort_by(|a, b| {
            b.1 .0
                .partial_cmp(&a.1 .0)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        for (strain_id, (abundance, confidence)) in strains {
            if *abundance > 1e-6 {
                report.push_str(&format!(
                    "  - {}: Abundance ~{:.2}%, Confidence ~{:.2}\n", // Use confidence proxy
                    strain_id,
                    abundance * 100.0,
                    confidence // Confidence is now 0-1 score
                ));
            }
        }
        report.push('\n');
    } else if results
        .classifications
        .first()
        .map_or(false, |c| c.level <= TaxonomicLevel::Species)
    {
        report.push_str(
            "Strain Abundance Estimates: No significant strain abundance detected or resolved.\n\n",
        );
    }

    // Footer
    report.push_str("----\n");
    report.push_str(&format!(
        "Results JSON: {}\n",
        results
            .results_file
            .as_ref()
            .map_or_else(|| "Not saved".to_string(), |p| p.display().to_string()) // Use map_or_else
    ));

    Ok(report)
}

// --- CLI Runner --- (Keeping definition as provided)

/// Command-line interface function to run the FASTQ processor.
pub fn run_fastq_cli(
    fastq_path: impl AsRef<Path>,
    sample_id: &str,
    db_path: impl AsRef<Path>,
    output_dir: impl AsRef<Path>,
    threads: usize,
) -> Result<(), ProcessingError> {
    // Initialize logger
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    // Log inputs
    info!("=== Starting FASTQ Processing ===");
    info!("Sample ID: {}", sample_id);
    info!("Input FASTQ: {}", fastq_path.as_ref().display());
    info!("Database Path: {}", db_path.as_ref().display());
    info!("Output Directory: {}", output_dir.as_ref().display());
    info!("Threads: {}", threads);

    // Hardcoded parameters (replace with CLI parsing in a real app)
    let macro_k = 31;
    let meso_k = 21;
    let sketch_size = 1000;
    let qc_params = QualityControlParams::default();
    info!(
        "Parameters: Macro K={}, Meso K={}, Sketch Size={}, QC={:?}",
        macro_k, meso_k, sketch_size, qc_params
    );

    // Create and initialize processor
    let mut processor = FastqProcessor::new(
        db_path,             // No need for as_ref() here, new() takes impl AsRef
        output_dir.as_ref(), // Use output_dir as cache_dir for simplicity
        threads,
        macro_k,
        meso_k,
        sketch_size,
        Some(qc_params),
        None, // No API key
    )?;

    processor.init_classifier()?; // Initialize before processing

    // Process the file
    let results = processor.process_file(
        fastq_path, // No need for as_ref() here
        sample_id, output_dir, // No need for as_ref() here
    )?;

    // Generate and print report
    info!("Generating final report...");
    let report = generate_report(&results)?;
    println!("\n{}", report); // Print report to stdout

    // Log output file location
    if let Some(ref report_file) = results.results_file {
        info!("Detailed results saved to: {}", report_file.display());
    } else {
        warn!("Results file path was not set."); // Should not happen if saved correctly
    }

    info!(
        "=== Processing Finished Successfully for Sample: {} ===",
        sample_id
    );
    Ok(())
}
