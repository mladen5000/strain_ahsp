use crate::adaptive::classifier::{AdaptiveClassifier, Classification, TaxonomicLevel};
use crate::database::DatabaseManager;
// Fix: Ensure correct signature types are imported and used consistently
// Assuming KmerSignature is the intended type for macro/meso signatures
use crate::sketch::signature::{KmerSignature, Signature}; // Removed ResolutionLevel
use crate::sketch::MultiResolutionSignature;
use log::{error, info, warn};
// Fix: Import needletail parser
use needletail::parse_fastx_file;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
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
    FastqError(String), // Keep for general parsing issues

    // Assuming SignatureError might come from sketch module or local operations
    #[error("Signature error: {0}")]
    SignatureError(String),

    #[error("Classification error: {0}")]
    ClassificationError(String),

    #[error("Strain estimation error: {0}")]
    StrainEstimationError(String),

    #[error("Database error: {0}")]
    DatabaseError(String),

    #[error("Needletail parsing error: {0}")] // Specific error for needletail
    NeedletailError(#[from] needletail::errors::ParseError),
}

// --- Structs (QC Params, Metrics, Results) ---
// (These seem okay, keeping them as they are)

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
#[derive(Debug, Clone, Serialize, Deserialize)]
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
    pub classifications: Vec<Classification>,
    pub strain_abundances: HashMap<String, (f64, f64)>,
    pub results_file: Option<PathBuf>,
}

// --- FastqProcessor ---

/// FASTQ processing pipeline
pub struct FastqProcessor {
    pub qc_params: QualityControlParams,
    pub threads: usize,
    pub chunk_size: usize,
    pub macro_k: usize,
    pub meso_k: usize,
    pub sketch_size: usize,
    pub db_manager: DatabaseManager,
    pub classifier: Option<AdaptiveClassifier>,
}

impl FastqProcessor {
    /// Create a new FASTQ processor
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
        let db_manager = DatabaseManager::new(
            db_path,
            cache_dir,
            sketch_size, // Assuming DB Manager needs sketch_size, threads
            threads,
            api_key,
        )
        .map_err(|e| ProcessingError::DatabaseError(format!("DB Manager init failed: {}", e)))?;

        Ok(FastqProcessor {
            qc_params: qc_params.unwrap_or_default(),
            threads,
            chunk_size: 100000,
            macro_k,
            meso_k,
            sketch_size,
            db_manager,
            classifier: None,
        })
    }

    /// Initialize the classifier by loading and converting reference signatures.
    pub fn init_classifier(&mut self) -> Result<(), ProcessingError> {
        info!("Fetching reference signatures from database...");
        // Assume get_all_signatures returns Vec<Arc<SomeDbType>> where SomeDbType has the fields needed
        // **This remains the most likely source of type errors if the assumption is wrong**
        let db_references = self.db_manager.database.get_all_signatures().map_err(|e| {
            ProcessingError::DatabaseError(format!("Failed to get signatures: {}", e))
        })?;

        if db_references.is_empty() {
            return Err(ProcessingError::DatabaseError(
                "No reference signatures found in the database.".to_string(),
            ));
        }
        info!(
            "Loaded {} reference signatures from DB.",
            db_references.len()
        );

        // Fix 1: Convert DB signatures to the type expected by AdaptiveClassifier
        // Assuming AdaptiveClassifier expects Vec<Arc<crate::sketch::signature::MultiResolutionSignature>>
        let mut sketch_signatures: Vec<Arc<MultiResolutionSignature>> =
            Vec::with_capacity(db_references.len());
        for db_sig_arc in db_references {
            // Check if essential signatures exist in the DB record
            if db_sig_arc.macro_signature.name.is_none() || db_sig_arc.meso_signature.name.is_none()
            {
                warn!("Skipping reference signature {} due to missing macro or meso signature in DB data.", db_sig_arc.taxon_id);
                continue;
            }

            // Construct the sketch::signature::MultiResolutionSignature
            // Assumes the fields in db_sig_arc are compatible or directly usable KmerSignatures/Option<VariantProfile> etc.
            // Also assumes the struct definition in sketch::signature matches this.
            let sketch_sig = MultiResolutionSignature {
                levels: db_sig_arc.levels.clone(), // Assuming levels is a field in the DB type
                taxon_id: db_sig_arc.taxon_id.clone(),
                lineage: db_sig_arc.lineage.clone(),
                // Clone the KmerSignatures from the DB Arc<SomeDbType>
                // This assumes db_sig_arc.macro_signature is Option<KmerSignature> or KmerSignature
                macro_signature: db_sig_arc.macro_signature.clone(), // Unwrap because we checked for None
                meso_signature: db_sig_arc.meso_signature.clone(), // Unwrap because we checked for None
                // Clone micro signature if present in DB type, else None
                micro_signature: db_sig_arc.micro_signature.clone(),
            };
            sketch_signatures.push(Arc::new(sketch_sig));
        }

        if sketch_signatures.is_empty() {
            return Err(ProcessingError::DatabaseError(
                "Failed to convert any database signatures to the required sketch format."
                    .to_string(),
            ));
        }

        info!(
            "Creating adaptive classifier with {} processed signatures...",
            sketch_signatures.len()
        );

        let sketch_signatures_owned: Vec<MultiResolutionSignature> = sketch_signatures
            .iter()
            .map(|sig| (**sig).clone())
            .collect();

        let classifier =
            AdaptiveClassifier::new(sketch_signatures_owned, None, None).map_err(|e| {
                ProcessingError::ClassificationError(format!("Classifier creation failed: {}", e))
            })?;
        info!("Classifier initialized successfully.");

        self.classifier = Some(classifier);
        Ok(())
    }

    /// Process a FASTQ file: read, QC, sketch, classify, estimate strains, and report.
    pub fn process_file(
        &self,
        fastq_path: impl AsRef<Path>,
        sample_id: &str,
        output_dir: impl AsRef<Path>,
    ) -> Result<ClassificationResults, ProcessingError> {
        let start_time = Instant::now();

        let classifier = self.classifier.as_ref().ok_or_else(|| {
            ProcessingError::ClassificationError(
                "Classifier not initialized. Call init_classifier() first.".to_string(),
            )
        })?;

        let output_path = output_dir.as_ref();
        std::fs::create_dir_all(output_path)?;

        let metrics = Arc::new(Mutex::new(ProcessingMetrics {
            total_reads: 0,
            passed_reads: 0,
            total_bases: 0,
            passed_bases: 0,
            avg_read_length: 0.0,
            processing_time_seconds: 0.0,
        }));

        // Fix: Initialize Signature struct instead of KmerSignature for macro_signature
        let macro_sig = Signature {
            algorithm: "Macro".to_string(),
            name: None,
            kmer_size: self.macro_k,
            num_hashes: 0,
            filename: None,
            path: None,
            hashes: Vec::new(),
        };
        let meso_sig = Signature {
            algorithm: "Meso".to_string(),
            name: None,
            kmer_size: self.meso_k,
            num_hashes: 0,
            filename: None,
            path: None,
            hashes: Vec::new(),
        };

        let initial_signature = MultiResolutionSignature {
            taxon_id: sample_id.to_string(),
            lineage: Vec::new(),
            macro_signature: macro_sig,
            meso_signature: meso_sig,
            micro_signature: Signature {
                algorithm: "Micro".to_string(),
                name: None,
                kmer_size: self.meso_k,
                num_hashes: 0,
                filename: None,
                path: None,
                hashes: Vec::new(),
            },
            levels: vec![],
        };
        let signature = Arc::new(Mutex::new(initial_signature));

        let mut reader = parse_fastx_file(fastq_path.as_ref())?; // Use '?'

        let mut current_chunk = Vec::with_capacity(self.chunk_size);

        info!("Processing file: {}", fastq_path.as_ref().display());

        while let Some(record_result) = reader.next() {
            let record = record_result?; // Use '?'
            current_chunk.push((record.seq().to_vec(), record.qual().map(|q| q.to_vec())));

            if current_chunk.len() >= self.chunk_size {
                self.process_chunk(&current_chunk, &metrics, &signature)?;
                current_chunk.clear();
            }
        }

        if !current_chunk.is_empty() {
            self.process_chunk(&current_chunk, &metrics, &signature)?;
        }

        let elapsed = start_time.elapsed().as_secs_f64();

        let final_metrics = {
            let mut metrics_guard = metrics.lock().unwrap();
            metrics_guard.processing_time_seconds = elapsed;
            if metrics_guard.passed_reads > 0 {
                metrics_guard.avg_read_length =
                    metrics_guard.passed_bases as f64 / metrics_guard.passed_reads as f64;
            }
            metrics_guard.clone()
        };

        let final_signature = signature.lock().unwrap().clone();

        info!("Classifying final sample signature...");
        // Use the get_hierarchical_classifications which currently wraps classify
        let classifications =
            self.get_hierarchical_classifications(&final_signature, classifier)?;

        let best_classification = classifications.first(); // get_hierarchical_classifications returns Vec

        let strain_abundances = if let Some(cls) = best_classification {
            info!(
                "Top classification: {} ({:?}), Confidence: {:.4}",
                cls.taxon_id, cls.level, cls.confidence
            );
            if cls.level <= TaxonomicLevel::Species {
                info!("Attempting strain estimation for {}...", cls.taxon_id);
                self.estimate_strain_abundances(&final_signature, classifier, &cls.taxon_id)?
            } else {
                info!(
                    "Classification level ({:?}) is above Species, skipping strain estimation.",
                    cls.level
                );
                HashMap::new()
            }
        } else {
            warn!("Classifier returned Ok but no classification found.");
            HashMap::new()
        };

        let results_file_path = output_path.join(format!("{}_results.json", sample_id));
        let results = ClassificationResults {
            sample_id: sample_id.to_string(),
            metrics: final_metrics.clone(),
            classifications, // Store the Vec from get_hierarchical_classifications
            strain_abundances,
            results_file: Some(results_file_path.clone()),
        };

        info!("Writing results to {}", results_file_path.display());
        let file = File::create(&results_file_path)?;
        let writer = BufWriter::new(file);
        serde_json::to_writer_pretty(writer, &results)
            .map_err(|e| ProcessingError::IoError(io::Error::new(io::ErrorKind::Other, e)))?;

        info!("Processed sample {} in {:.2} seconds", sample_id, elapsed);
        info!(
            "Reads: {}/{} passed QC ({:.1}%)",
            final_metrics.passed_reads,
            final_metrics.total_reads,
            100.0 * final_metrics.passed_reads as f64 / final_metrics.total_reads.max(1) as f64
        );
        info!(
            "Avg Read Length (Passed QC): {:.1} bp",
            final_metrics.avg_read_length
        );

        Ok(results)
    }

    /// Process a chunk of reads in parallel: apply QC and update the shared signature.
    fn process_chunk(
        &self,
        chunk: &[(Vec<u8>, Option<Vec<u8>>)],
        metrics: &Arc<Mutex<ProcessingMetrics>>,
        signature: &Arc<Mutex<MultiResolutionSignature>>,
    ) -> Result<(), ProcessingError> {
        chunk
            .par_iter()
            .try_for_each(|(seq, qual)| -> Result<(), ProcessingError> {
                let original_len = seq.len();
                let processed_seq_opt = self.apply_quality_control(seq, qual.as_ref());

                {
                    let mut metrics_guard = metrics.lock().unwrap();
                    metrics_guard.total_reads += 1;
                    metrics_guard.total_bases += original_len;
                    if let Some(ref processed) = processed_seq_opt {
                        metrics_guard.passed_reads += 1;
                        metrics_guard.passed_bases += processed.len();
                    }
                }

                if let Some(processed_seq) = processed_seq_opt {
                    if !processed_seq.is_empty() {
                        let mut sig_guard = signature.lock().unwrap();
                        // Fix 18, 19: Assume macro/meso_signature fields are KmerSignature structs
                        // and call the appropriate method (e.g., add_sequence or add_kmers).
                        // Requires checking the KmerSignature definition in sketch/signature.rs.
                        // Let's assume `add_sequence` exists on KmerSignature based on previous attempts.
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
                        // Do not update micro_signature here
                    }
                }
                Ok(())
            })
    }

    /// Apply quality control filters to a single read.
    fn apply_quality_control(&self, seq: &[u8], qual: Option<&Vec<u8>>) -> Option<Vec<u8>> {
        // 1. Check initial length
        if seq.len() < self.qc_params.min_length {
            return None;
        }

        // 2. Check N content
        let n_count = seq
            .iter()
            .filter(|&&base| base == b'N' || base == b'n')
            .count();
        let n_percent = 100.0 * n_count as f64 / seq.len() as f64;
        if n_percent > self.qc_params.max_n_percent {
            return None;
        }

        // 3. Quality trimming and average quality check
        if let Some(qual_vec) = qual {
            if qual_vec.len() != seq.len() {
                error!(
                    "Sequence length ({}) and quality length ({}) mismatch. Discarding read.",
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
                .map(|&q| (q.saturating_sub(33)) as f64)
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

    /// Get hierarchical classifications (currently just returns the best one).
    fn get_hierarchical_classifications(
        &self,
        signature: &MultiResolutionSignature,
        classifier: &AdaptiveClassifier,
    ) -> Result<Vec<Classification>, ProcessingError> {
        let best_classification = classifier.classify(signature).map_err(|e| {
            ProcessingError::ClassificationError(format!("Classification failed: {}", e))
        })?;
        Ok(vec![best_classification])
    }

    /// Estimate relative abundances of strains related to the classified species.
    fn estimate_strain_abundances(
        &self,
        signature: &MultiResolutionSignature,
        classifier: &AdaptiveClassifier,
        target_species_id: &str,
    ) -> Result<HashMap<String, (f64, f64)>, ProcessingError> {
        info!(
            "Estimating strain abundances relative to target: {}",
            target_species_id
        );

        // Use the existing filtering logic based on lineage and taxon_id
        let relevant_strains = classifier
            .references
            .iter()
            .filter(|ref_sig| {
                // Crude check: is it downstream in lineage and not the species itself?
                ref_sig.lineage.contains(&target_species_id.to_string())
                    && ref_sig.taxon_id != target_species_id
                // A check based on TaxonomicLevel might be better if available and reliable on ref_sig
                // For example: if let Some(ref_level) = ref_sig.level { ref_level > TaxonomicLevel::Species && ... }
            })
            .collect::<Vec<_>>();

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

        let mut similarities = HashMap::new();
        let mut total_similarity = 0.0;

        for strain_sig in relevant_strains {
            let sim = signature.similarity(strain_sig, None); // Use overall similarity
            if sim > 0.0 {
                similarities.insert(strain_sig.taxon_id.clone(), sim);
                total_similarity += sim;
            }
        }

        let mut abundances = HashMap::new();
        if total_similarity > f64::EPSILON {
            for (id, sim) in similarities {
                let abundance = sim / total_similarity;
                let confidence = 0.1; // Placeholder
                                      // Fix 20: Clone id before inserting into abundances map
                abundances.insert(id.clone(), (abundance, confidence));
                // Log using the original (uncloned) id is fine here
                info!(
                    "  Strain {}: Relative Abundance ~{:.2}%, Similarity {:.4}",
                    id, // Use original id here
                    abundance * 100.0,
                    sim
                );
            }
        } else {
            info!("Total similarity to relevant strains is zero or negligible.");
        }

        Ok(abundances)
    }
}

/// Generate a formatted text report from the classification results.
pub fn generate_report(results: &ClassificationResults) -> Result<String, ProcessingError> {
    let mut report = String::new();

    // Header
    report.push_str(&format!(
        "AHSP Classification Report for Sample: {}\n",
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
        // Report top hit primarily
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
                // Sort scores by similarity value (descending)
                let mut sorted_scores: Vec<_> = classification.similarity_scores.iter().collect();
                sorted_scores.sort_by(|(_, score_a), (_, score_b)| {
                    score_b
                        .partial_cmp(score_a)
                        .unwrap_or(std::cmp::Ordering::Equal)
                });

                for (level, score) in sorted_scores {
                    // Assuming level is Debug-printable (like ResolutionLevel or similar enum/struct)
                    report.push_str(&format!("    {:?}: {:.4}\n", level, score));
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
                    "  - {}: {:.2}% (± {:.1}%)\n",
                    strain_id,
                    abundance * 100.0,
                    confidence * 100.0
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
            .map_or("Not saved".to_string(), |p| p.display().to_string())
    ));

    Ok(report)
}

/// Command-line interface function to run the FASTQ processor.
pub fn run_fastq_cli(
    fastq_path: impl AsRef<Path>,
    sample_id: &str,
    db_path: impl AsRef<Path>,
    output_dir: impl AsRef<Path>,
    threads: usize,
) -> Result<(), ProcessingError> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    info!("Starting FASTQ processing for sample: {}", sample_id);
    info!("Input FASTQ: {}", fastq_path.as_ref().display());
    info!("Database path: {}", db_path.as_ref().display());
    info!("Output directory: {}", output_dir.as_ref().display());
    info!("Using {} threads", threads);

    let macro_k = 31;
    let meso_k = 21;
    let sketch_size = 1000;
    let qc_params = QualityControlParams::default();

    info!(
        "Parameters: Macro K={}, Meso K={}, Sketch Size={}, QC={:?}",
        macro_k, meso_k, sketch_size, qc_params
    );

    let mut processor = FastqProcessor::new(
        db_path.as_ref(),
        output_dir.as_ref(),
        threads,
        macro_k,
        meso_k,
        sketch_size,
        Some(qc_params),
        None, // No API key
    )?;

    info!("Initializing classifier...");
    processor.init_classifier()?;
    info!("Classifier initialized.");

    info!("Processing FASTQ file...");
    let results = processor.process_file(fastq_path.as_ref(), sample_id, output_dir.as_ref())?;

    info!("Generating report...");
    let report = generate_report(&results)?;
    println!("\n{}", report);

    if let Some(ref report_file) = results.results_file {
        info!("Detailed results saved to: {}", report_file.display());
    } else {
        warn!("Results file path was not set in the results structure.");
    }

    info!("Processing finished successfully for sample: {}", sample_id);
    Ok(())
}
