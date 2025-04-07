use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::time::Instant;

use log::{error, info, warn};
use needletail::{parse_fastx_file, Sequence};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use thiserror::Error;

use crate::adaptive::classifier::{AdaptiveClassifier, Classification, TaxonomicLevel};
use crate::bio::signature::{
    KmerSignature, MultiResolutionSignature, ResolutionLevel, VariantProfile,
};
use crate::database::{DatabaseManager, SignatureDatabase};
use crate::stats::bayesian::StrainMixtureModel;

#[derive(Error, Debug)]
pub enum ProcessingError {
    #[error("IO error: {0}")]
    IoError(#[from] io::Error),

    #[error("FASTQ parsing error: {0}")]
    FastqError(String),

    #[error("Signature error: {0}")]
    SignatureError(String),

    #[error("Classification error: {0}")]
    ClassificationError(String),

    #[error("Strain estimation error: {0}")]
    StrainEstimationError(String),

    #[error("Database error: {0}")]
    DatabaseError(String),
}

/// Quality control parameters for FASTQ processing
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityControlParams {
    /// Minimum average quality score
    pub min_avg_quality: f64,

    /// Minimum read length after trimming
    pub min_length: usize,

    /// Quality score threshold for trimming
    pub trim_quality: u8,

    /// Maximum percentage of N bases allowed
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
    /// Total number of reads processed
    pub total_reads: usize,

    /// Number of reads passing QC
    pub passed_reads: usize,

    /// Total bases processed
    pub total_bases: usize,

    /// Number of bases passing QC
    pub passed_bases: usize,

    /// Average read length
    pub avg_read_length: f64,

    /// Processing time in seconds
    pub processing_time_seconds: f64,
}

/// Sample classification results
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ClassificationResults {
    /// Sample ID
    pub sample_id: String,

    /// Processing metrics
    pub metrics: ProcessingMetrics,

    /// Top taxonomic classifications
    pub classifications: Vec<Classification>,

    /// Strain abundance estimates
    pub strain_abundances: HashMap<String, (f64, f64)>, // strain_id -> (abundance, confidence)

    /// Path to results file
    pub results_file: Option<PathBuf>,
}

/// FASTQ processing pipeline
pub struct FastqProcessor {
    /// Quality control parameters
    pub qc_params: QualityControlParams,

    /// Number of threads to use
    pub threads: usize,

    /// Chunk size for parallel processing
    pub chunk_size: usize,

    /// K-mer size for macro-level signatures
    pub macro_k: usize,

    /// K-mer size for meso-level signatures
    pub meso_k: usize,

    /// Sketch size for signatures
    pub sketch_size: usize,

    /// Database manager
    pub db_manager: DatabaseManager,

    /// Classifier
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
        // Create database manager
        let db_manager = DatabaseManager::new(
            db_path,
            cache_dir,
            macro_k,
            meso_k,
            sketch_size,
            threads,
            api_key,
        )
        .map_err(|e| ProcessingError::DatabaseError(format!("{}", e)))?;

        Ok(FastqProcessor {
            qc_params: qc_params.unwrap_or_default(),
            threads,
            chunk_size: 100000, // Process 100k reads at a time
            macro_k,
            meso_k,
            sketch_size,
            db_manager,
            classifier: None, // Will be initialized later
        })
    }

    /// Initialize the classifier
    pub fn init_classifier(&mut self) -> Result<(), ProcessingError> {
        // Get all reference signatures
        let references = self
            .db_manager
            .database
            .get_all_signatures()
            .map_err(|e| ProcessingError::DatabaseError(format!("{}", e)))?;

        // Create adaptive classifier
        let classifier = AdaptiveClassifier::new(references, None, None)
            .map_err(|e| ProcessingError::ClassificationError(format!("{}", e)))?;

        self.classifier = Some(classifier);

        Ok(())
    }

    /// Process a FASTQ file
    pub fn process_file(
        &self,
        fastq_path: impl AsRef<Path>,
        sample_id: &str,
        output_dir: impl AsRef<Path>,
    ) -> Result<ClassificationResults, ProcessingError> {
        let start_time = Instant::now();

        // Check if classifier is initialized
        let classifier = match &self.classifier {
            Some(c) => c,
            None => {
                return Err(ProcessingError::ClassificationError(
                    "Classifier not initialized. Call init_classifier() first.".to_string(),
                ))
            }
        };

        // Create output directory if it doesn't exist
        let output_path = output_dir.as_ref();
        if !output_path.exists() {
            std::fs::create_dir_all(output_path)?;
        }

        // Initialize metrics
        let metrics = Arc::new(Mutex::new(ProcessingMetrics {
            total_reads: 0,
            passed_reads: 0,
            total_bases: 0,
            passed_bases: 0,
            avg_read_length: 0.0,
            processing_time_seconds: 0.0,
        }));

        // Initialize signature
        let signature = Arc::new(Mutex::new(
            MultiResolutionSignature::new(
                sample_id,
                Vec::new(),
                self.macro_k,
                self.meso_k,
                self.sketch_size,
                self.sketch_size,
            )
            .map_err(|e| ProcessingError::SignatureError(format!("{}", e)))?,
        ));

        // Create a reader for the FASTQ file
        let mut reader = parse_fastx_file(fastq_path.as_ref())
            .map_err(|e| ProcessingError::FastqError(format!("{}", e)))?;

        // Read and process chunks of reads in parallel
        let mut current_chunk = Vec::new();

        info!("Processing file: {}", fastq_path.as_ref().display());

        while let Some(record) = reader.next() {
            let record = record.map_err(|e| ProcessingError::FastqError(format!("{}", e)))?;

            current_chunk.push((record.seq().to_vec(), record.qual().map(|q| q.to_vec())));

            // Process chunk when it reaches the specified size
            if current_chunk.len() >= self.chunk_size {
                self.process_chunk(&current_chunk, &metrics, &signature)?;
                current_chunk.clear();
            }
        }

        // Process any remaining reads
        if !current_chunk.is_empty() {
            self.process_chunk(&current_chunk, &metrics, &signature)?;
        }

        // Calculate processing time
        let elapsed = start_time.elapsed().as_secs_f64();

        // Update metrics
        let mut metrics = metrics.lock().unwrap();
        metrics.processing_time_seconds = elapsed;

        if metrics.passed_reads > 0 {
            metrics.avg_read_length = metrics.passed_bases as f64 / metrics.passed_reads as f64;
        }

        // Get the signature
        let signature = signature.lock().unwrap().clone();

        // Classify the sample
        let classification = classifier
            .classify(&signature)
            .map_err(|e| ProcessingError::ClassificationError(format!("{}", e)))?;

        // Get top-level classifications at different taxonomic levels
        let classifications = self.get_hierarchical_classifications(&signature, classifier)?;

        // Estimate strain abundances if we have species-level classification
        let strain_abundances = if classification.level <= TaxonomicLevel::Species {
            self.estimate_strain_abundances(&signature, classifier)?
        } else {
            HashMap::new()
        };

        // Write results to file
        let results_file = output_path.join(format!("{}_results.json", sample_id));
        let results = ClassificationResults {
            sample_id: sample_id.to_string(),
            metrics: metrics.clone(),
            classifications,
            strain_abundances,
            results_file: Some(results_file.clone()),
        };

        // Serialize and write results
        let file = File::create(&results_file)?;
        let writer = BufWriter::new(file);
        serde_json::to_writer_pretty(writer, &results)
            .map_err(|e| ProcessingError::IoError(io::Error::new(io::ErrorKind::Other, e)))?;

        info!("Processed sample {} in {:.2} seconds", sample_id, elapsed);
        info!(
            "Reads: {}/{} passed QC ({:.1}%)",
            metrics.passed_reads,
            metrics.total_reads,
            100.0 * metrics.passed_reads as f64 / metrics.total_reads.max(1) as f64
        );

        Ok(results)
    }

    /// Process a chunk of reads in parallel
    fn process_chunk(
        &self,
        chunk: &[(Vec<u8>, Option<Vec<u8>>)],
        metrics: &Arc<Mutex<ProcessingMetrics>>,
        signature: &Arc<Mutex<MultiResolutionSignature>>,
    ) -> Result<(), ProcessingError> {
        // Process reads in parallel
        let results: Vec<Option<Vec<u8>>> = chunk
            .par_iter()
            .map(|(seq, qual)| {
                // Apply quality control
                if let Some(processed) = self.apply_quality_control(seq, qual.as_ref()) {
                    // Update metrics
                    let mut metrics_guard = metrics.lock().unwrap();
                    metrics_guard.total_reads += 1;
                    metrics_guard.total_bases += seq.len();
                    metrics_guard.passed_reads += 1;
                    metrics_guard.passed_bases += processed.len();

                    Some(processed)
                } else {
                    // Read failed QC
                    let mut metrics_guard = metrics.lock().unwrap();
                    metrics_guard.total_reads += 1;
                    metrics_guard.total_bases += seq.len();

                    None
                }
            })
            .collect();

        // Add passing reads to signature
        let mut sig_guard = signature.lock().unwrap();
        for seq_opt in results {
            if let Some(seq) = seq_opt {
                sig_guard
                    .add_sequence(&seq)
                    .map_err(|e| ProcessingError::SignatureError(format!("{}", e)))?;
            }
        }

        Ok(())
    }

    /// Apply quality control to a read
    fn apply_quality_control(&self, seq: &[u8], qual: Option<&Vec<u8>>) -> Option<Vec<u8>> {
        // Check length
        if seq.len() < self.qc_params.min_length {
            return None;
        }

        // Check N content
        let n_count = seq.iter().filter(|&&base| base == b'N').count();
        let n_percent = 100.0 * n_count as f64 / seq.len() as f64;
        if n_percent > self.qc_params.max_n_percent {
            return None;
        }

        // Check and apply quality trimming
        if let Some(qual_vec) = qual {
            // Calculate average quality
            let avg_quality =
                qual_vec.iter().map(|&q| (q - 33) as f64).sum::<f64>() / qual_vec.len() as f64;

            if avg_quality < self.qc_params.min_avg_quality {
                return None;
            }

            // Find quality-based trim positions
            let mut trim_start = 0;
            let mut trim_end = seq.len();

            // Trim from start
            for (i, &q) in qual_vec.iter().enumerate() {
                if (q - 33) as u8 >= self.qc_params.trim_quality {
                    trim_start = i;
                    break;
                }
            }

            // Trim from end
            for i in (0..qual_vec.len()).rev() {
                if (qual_vec[i] - 33) as u8 >= self.qc_params.trim_quality {
                    trim_end = i + 1;
                    break;
                }
            }

            // Check if trimmed read is long enough
            if trim_end - trim_start < self.qc_params.min_length {
                return None;
            }

            // Return trimmed sequence
            Some(seq[trim_start..trim_end].to_vec())
        } else {
            // No quality scores available, return the original sequence
            Some(seq.to_vec())
        }
    }

    /// Get hierarchical classifications at different taxonomic levels
    fn get_hierarchical_classifications(
        &self,
        signature: &MultiResolutionSignature,
        classifier: &AdaptiveClassifier,
    ) -> Result<Vec<Classification>, ProcessingError> {
        // Get classification at different resolution levels
        let macro_classification = classifier
            .classify(signature)
            .map_err(|e| ProcessingError::ClassificationError(format!("{}", e)))?;

        // Create a vector of classifications at different levels
        let mut classifications = Vec::new();
        classifications.push(macro_classification);

        // Add any additional hierarchical classifications if needed
        // For example, if we classified at species level, we might want to add genus, family, etc.

        Ok(classifications)
    }

    /// Estimate strain abundances in the sample
    fn estimate_strain_abundances(
        &self,
        signature: &MultiResolutionSignature,
        classifier: &AdaptiveClassifier,
    ) -> Result<HashMap<String, (f64, f64)>, ProcessingError> {
        // For simplicity, we'll use the best-matching strains from the classification
        // In a real implementation, we would extract k-mer profiles and use the Bayesian model

        // Get related reference signatures
        let species_id = signature.taxon_id.clone();
        let related_strains = classifier
            .references
            .iter()
            .filter(|ref_sig| {
                ref_sig.lineage.contains(&species_id) || ref_sig.taxon_id == species_id
            })
            .collect::<Vec<_>>();

        if related_strains.is_empty() {
            return Ok(HashMap::new());
        }

        // Calculate similarities to each related strain
        let mut similarities = HashMap::new();
        let mut total_similarity = 0.0;

        for strain in related_strains {
            let sim = signature.similarity(strain, None);
            similarities.insert(strain.taxon_id.clone(), sim);
            total_similarity += sim;
        }

        // Normalize to get abundance estimates
        let mut abundances = HashMap::new();

        if total_similarity > 0.0 {
            for (id, sim) in similarities {
                // Normalize and add confidence interval (simplified)
                let abundance = sim / total_similarity;
                let confidence = 0.1 * abundance; // Simplified confidence calculation
                abundances.insert(id, (abundance, confidence));
            }
        }

        Ok(abundances)
    }
}

/// Generate a text report from classification results
pub fn generate_report(results: &ClassificationResults) -> Result<String, ProcessingError> {
    let mut report = String::new();

    // Add header
    report.push_str(&format!(
        "AHSP Classification Report for Sample: {}\n",
        results.sample_id
    ));
    report.push_str("=================================================\n\n");

    // Add processing metrics
    report.push_str("Processing Metrics:\n");
    report.push_str(&format!("  Total reads: {}\n", results.metrics.total_reads));
    report.push_str(&format!(
        "  Passed QC: {} ({:.1}%)\n",
        results.metrics.passed_reads,
        100.0 * results.metrics.passed_reads as f64 / results.metrics.total_reads.max(1) as f64
    ));
    report.push_str(&format!(
        "  Average read length: {:.1} bp\n",
        results.metrics.avg_read_length
    ));
    report.push_str(&format!(
        "  Processing time: {:.2} seconds\n\n",
        results.metrics.processing_time_seconds
    ));

    // Add classification results
    report.push_str("Classification Results:\n");

    for (i, classification) in results.classifications.iter().enumerate() {
        report.push_str(&format!("Classification #{}\n", i + 1));
        report.push_str(&format!("  Taxon ID: {}\n", classification.taxon_id));
        report.push_str(&format!("  Taxonomic level: {:?}\n", classification.level));
        report.push_str(&format!("  Confidence: {:.2}\n", classification.confidence));

        // Add lineage
        if !classification.lineage.is_empty() {
            report.push_str("  Lineage: ");
            for (j, taxon) in classification.lineage.iter().enumerate() {
                if j > 0 {
                    report.push_str(" > ");
                }
                report.push_str(taxon);
            }
            report.push('\n');
        }

        // Add similarity scores
        report.push_str("  Similarity scores:\n");
        for (level, score) in &classification.similarity_scores {
            report.push_str(&format!("    {:?}: {:.2}\n", level, score));
        }
        report.push('\n');
    }

    // Add strain abundances if available
    if !results.strain_abundances.is_empty() {
        report.push_str("Strain Abundances:\n");

        // Sort strains by abundance
        let mut strains: Vec<_> = results.strain_abundances.iter().collect();
        strains.sort_by(|a, b| b.1 .0.partial_cmp(&a.1 .0).unwrap());

        for (strain_id, (abundance, confidence)) in strains {
            report.push_str(&format!(
                "  {}: {:.2}% (±{:.2}%)\n",
                strain_id,
                abundance * 100.0,
                confidence * 100.0
            ));
        }
        report.push('\n');
    }

    // Add footer
    report.push_str(&format!(
        "Results file: {}\n",
        results
            .results_file
            .as_ref()
            .map_or("None".to_string(), |p| p.display().to_string())
    ));

    Ok(report)
}

/// Command-line interface for FASTQ processing
pub fn run_fastq_cli(
    fastq_path: impl AsRef<Path>,
    sample_id: &str,
    db_path: impl AsRef<Path>,
    output_dir: impl AsRef<Path>,
    threads: usize,
) -> Result<(), ProcessingError> {
    // Initialize logger
    env_logger::init();

    // Create FASTQ processor
    let mut processor = FastqProcessor::new(
        db_path,
        "genome_cache", // Default cache directory
        threads,
        31,   // Default macro_k
        21,   // Default meso_k
        1000, // Default sketch_size
        None, // Default QC parameters
        None, // No API key
    )?;

    // Initialize classifier
    processor.init_classifier()?;

    // Process FASTQ file
    let results = processor.process_file(fastq_path, sample_id, output_dir)?;

    // Generate and print report
    let report = generate_report(&results)?;
    println!("{}", report);

    Ok(())
}
