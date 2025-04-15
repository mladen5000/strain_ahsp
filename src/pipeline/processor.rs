// #![allow(dead_code)] // Allow unused code during development/refinement

// // --- Existing Imports ---
// use crate::adaptive::classifier::Classification;
// use crate::bio::taxonomy::TaxonomicLineage;
// use crate::sketch::signature::{KmerSignature, MultiResolutionSignature, Signature};
// use log::{info, warn};
// use needletail::Sequence;
// use rayon::prelude::*;
// use serde::{Deserialize, Serialize};
// use std::collections::hash_map::DefaultHasher;
// use std::sync::{Arc, Mutex};
// use std::{
//     collections::HashMap,
//     fs,
//     hash::{Hash, Hasher},
//     io,
//     path::{Path, PathBuf},
//     time::Instant,
// };

// // --- Error Type ---
// #[derive(Debug)]
// pub enum ProcessingError {
//     IoError(io::Error),
//     ParseError(needletail::errors::ParseError),
//     ProcessError(String),
//     SignatureError(String),
// }

// impl From<io::Error> for ProcessingError {
//     fn from(err: io::Error) -> Self {
//         ProcessingError::IoError(err)
//     }
// }

// impl From<needletail::errors::ParseError> for ProcessingError {
//     fn from(err: needletail::errors::ParseError) -> Self {
//         ProcessingError::ParseError(err)
//     }
// }

// // --- Sketch Trait ---
// pub trait Sketch {
//     fn add_sequence(&mut self, seq: &[u8]) -> Result<(), String>;
//     fn len(&self) -> usize;
//     fn is_empty(&self) -> bool;
//     fn sketch_size(&self) -> usize;
//     fn best_level_similarity(&self, other: &Self) -> Result<f64, String>;
// }

// impl Sketch for MultiResolutionSignature {
//     fn add_sequence(&mut self, seq: &[u8]) -> Result<(), String> {
//         // Add sequence to all resolution levels
//         for level in &mut self.levels {
//             level.add_sequence(seq)?;
//         }
//         Ok(())
//     }

//     fn len(&self) -> usize {
//         // Use first level (macro) for length
//         self.levels.first().map_or(0, |s| s.sketch.hashes.len())
//     }

//     fn is_empty(&self) -> bool {
//         self.levels.iter().all(|s| s.sketch.hashes.is_empty())
//     }

//     fn sketch_size(&self) -> usize {
//         // Use first level for sketch size
//         self.levels.first().map_or(0, |s| s.sketch.num_hashes)
//     }

//     fn best_level_similarity(&self, other: &Self) -> Result<f64, String> {
//         // Calculate similarity at each level and return the maximum
//         let mut max_sim = 0.0;
//         for (self_level, other_level) in self.levels.iter().zip(other.levels.iter()) {
//             if let Some(sim) = self_level.jaccard_similarity(other_level) {
//                 max_sim = f64::max(max_sim, sim);
//             }
//         }
//         if max_sim.is_nan() {
//             return Err("Similarity calculation resulted in NaN".to_string());
//         }
//         Ok(max_sim)
//     }
// }

// impl Sketch for Signature {
//     fn add_sequence(&mut self, seq: &[u8]) -> Result<(), String> {
//         let k = self.k();
//         if seq.len() < k {
//             return Err("Sequence length is smaller than k-mer size".to_string());
//         }

//         let mut hasher = DefaultHasher::new();
//         for i in 0..=seq.len() - k {
//             let kmer = &seq[i..i + k];
//             ByteSlice(kmer).hash(&mut hasher);
//             let hash = hasher.finish();

//             if self.hashes.len() < self.num_hashes {
//                 self.hashes.push(hash);
//                 if self.hashes.len() == self.num_hashes {
//                     self.hashes.sort_unstable();
//                 }
//             } else if hash < self.hashes[self.num_hashes - 1] {
//                 self.hashes[self.num_hashes - 1] = hash;
//                 self.hashes.sort_unstable();
//             }
//         }
//         Ok(())
//     }

//     fn len(&self) -> usize {
//         self.hashes.len()
//     }

//     fn is_empty(&self) -> bool {
//         self.hashes.is_empty()
//     }

//     fn sketch_size(&self) -> usize {
//         self.num_hashes
//     }

//     fn best_level_similarity(&self, other: &Self) -> Result<f64, String> {
//         if self.is_empty() || other.is_empty() {
//             return Ok(0.0);
//         }

//         let intersection: Vec<_> = self
//             .hashes
//             .iter()
//             .filter(|h| other.hashes.contains(h))
//             .collect();

//         Ok(intersection.len() as f64 / self.hashes.len().max(other.hashes.len()) as f64)
//     }
// }

// // --- Implementing required traits ---
// impl Sketch for KmerSignature {
//     fn add_sequence(&mut self, seq: &[u8]) -> Result<(), String> {
//         // Process k-mers from sequence and add to kmers
//         let mut hasher = DefaultHasher::new();
//         for i in 0..=seq.len().saturating_sub(self.kmer_size) {
//             let kmer = &seq[i..i + self.kmer_size];
//             ByteSlice(kmer).hash(&mut hasher);
//             let hash = hasher.finish();
//             self.kmers.push(hash);
//         }
//         Ok(())
//     }

//     fn len(&self) -> usize {
//         self.kmers.len()
//     }

//     fn is_empty(&self) -> bool {
//         self.kmers.is_empty()
//     }

//     fn sketch_size(&self) -> usize {
//         self.num_hashes
//     }

//     fn best_level_similarity(&self, other: &Self) -> Result<f64, String> {
//         let intersection: Vec<_> = self
//             .kmers
//             .iter()
//             .filter(|h| other.kmers.contains(h))
//             .collect();

//         if self.kmers.is_empty() || other.kmers.is_empty() {
//             return Ok(0.0);
//         }

//         Ok(intersection.len() as f64 / self.kmers.len().max(other.kmers.len()) as f64)
//     }
// }

// impl From<KmerSignature> for Signature {
//     fn from(kmer_sig: KmerSignature) -> Self {
//         Self {
//             algorithm: "minhash".to_string(),
//             kmer_size: kmer_sig.kmer_size,
//             num_hashes: kmer_sig.num_hashes,
//             name: None,
//             filename: None,
//             path: None,
//             hashes: kmer_sig.kmers,
//         }
//     }
// }

// impl TryFrom<Signature> for KmerSignature {
//     type Error = String;

//     fn try_from(sig: Signature) -> Result<Self, Self::Error> {
//         Ok(Self {
//             kmers: sig.hashes,
//             counts: Vec::new(),
//             total_kmers: 0,
//             kmer_size: sig.kmer_size,
//             num_hashes: sig.num_hashes,
//         })
//     }
// }

// // Conversion from Mutex PoisonError
// impl<T> From<std::sync::PoisonError<T>> for ProcessingError {
//     fn from(err: std::sync::PoisonError<T>) -> Self {
//         ProcessingError::ProcessError(format!("Mutex poisoned: {}", err))
//     }
// }

// // --- Structs (QC Params, Metrics, Results) --- (No changes)

// #[derive(Debug, Clone, Serialize, Deserialize)]
// pub struct QualityControlParams {
//     pub min_avg_quality: f64,
//     pub min_length: usize,
//     pub trim_quality: u8,
//     pub max_n_percent: f64,
// }

// impl Default for QualityControlParams {
//     fn default() -> Self {
//         Self {
//             min_avg_quality: 20.0,
//             min_length: 50,
//             trim_quality: 20,
//             max_n_percent: 10.0,
//         }
//     }
// }

// #[derive(Debug, Clone, Serialize, Deserialize)]
// pub struct ProcessingMetrics {
//     pub total_reads: u64,
//     pub reads_passed_qc: u64,
//     pub avg_read_length: f64,
//     pub processing_time_seconds: f64,
//     pub strain_abundances: HashMap<String, (f64, f64)>,
//     pub processed_reads: u64,
//     pub processed_bases: usize,
// }

// impl Default for ProcessingMetrics {
//     fn default() -> Self {
//         Self {
//             total_reads: 0,
//             reads_passed_qc: 0,
//             avg_read_length: 0.0,
//             processing_time_seconds: 0.0,
//             strain_abundances: HashMap::new(),
//             processed_reads: 0,
//             processed_bases: 0,
//         }
//     }
// }

// #[derive(Debug, Clone, Serialize, Deserialize)]
// pub struct SketchSizeMetrics {
//     pub macro_sketch_size: usize,
//     pub meso_sketch_size: usize,
//     pub micro_sketch_size: Option<usize>,
// }

// #[derive(Debug, Clone, Serialize, Deserialize)]
// pub struct ClassificationResults {
//     pub sample_id: String,
//     pub metrics: ProcessingMetrics,
//     pub classifications: Vec<Classification>,
//     pub strain_abundances: HashMap<String, (f64, f64)>,
//     pub results_file: Option<PathBuf>,
// }

// #[derive(Debug, Serialize, Deserialize)]
// pub struct ProcessingResult {
//     pub metrics: ProcessingMetrics,
//     pub signatures: Vec<Arc<MultiResolutionSignature>>,
//     pub sample_id: String,
// }

// // --- Processor ---
// pub struct Processor;

// impl Processor {
//     pub fn new() -> Self {
//         Self {}
//     }

//     fn process_signatures(
//         &mut self,
//         signatures: Vec<MultiResolutionSignature>,
//     ) -> Result<Vec<Arc<MultiResolutionSignature>>, ProcessingError> {
//         let processed: Vec<Arc<MultiResolutionSignature>> =
//             signatures.into_iter().map(|sig| Arc::new(sig)).collect();
//         Ok(processed)
//     }

//     fn create_empty_sample_signature(
//         &self,
//         sample_id: &str,
//     ) -> Result<Arc<Mutex<MultiResolutionSignature>>, ProcessingError> {
//         let mut signature = MultiResolutionSignature::new(sample_id.to_string(), Vec::new());

//         // Create three resolution levels with decreasing k-mer sizes and increasing sensitivity
//         let macro_k = 21;
//         let meso_k = 15;
//         let micro_k = 11;

//         let macro_sig = KmerSignature::new(
//             macro_k,
//             "DNA".to_string(),
//             "minhash".to_string(),
//             1000, // Larger sketch for coarse resolution
//             0,
//         );
//         let meso_sig = KmerSignature::new(
//             meso_k,
//             "DNA".to_string(),
//             "minhash".to_string(),
//             2000, // Medium sketch size
//             0,
//         );
//         let micro_sig = KmerSignature::new(
//             micro_k,
//             "DNA".to_string(),
//             "minhash".to_string(),
//             4000, // Smaller sketch for fine resolution
//             0,
//         );

//         signature.add_level(macro_sig);
//         signature.add_level(meso_sig);
//         signature.add_level(micro_sig);

//         Ok(Arc::new(Mutex::new(signature)))
//     }

//     fn log_summary_metrics(&self, metrics: &ProcessingMetrics) {
//         info!("Processing summary:");
//         info!("Total reads processed: {}", metrics.total_reads);
//         info!("Reads passed QC: {}", metrics.reads_passed_qc);
//         info!("Average read length: {:.2}", metrics.avg_read_length);
//     }

//     pub fn process_file(
//         &self,
//         fastq_path: impl AsRef<Path>,
//         sample_id: String,
//         output_dir: impl AsRef<Path>,
//     ) -> Result<ProcessingResult, ProcessingError> {
//         let start_time = Instant::now();
//         let mut metrics = ProcessingMetrics::default();
//         let mut signatures = Vec::new();

//         let mut reader = needletail::parse_fastx_file(fastq_path.as_ref())?;
//         while let Some(record) = reader.next() {
//             let record = record?;
//             metrics.total_reads += 1;

//             if let Some(processed_seq) = self.apply_quality_control(record.sequence(), None) {
//                 metrics.reads_passed_qc += 1;
//                 metrics.avg_read_length += processed_seq.len() as f64;
//                 // Process sequence...
//             }
//         }

//         if metrics.reads_passed_qc > 0 {
//             metrics.avg_read_length /= metrics.reads_passed_qc as f64;
//         }
//         metrics.processing_time_seconds = start_time.elapsed().as_secs_f64();

//         let processed_sigs: Vec<Arc<MultiResolutionSignature>> = signatures
//             .into_iter()
//             .map(|sig| {
//                 Arc::new(MultiResolutionSignature {
//                     taxon_id: sig.taxon_id.unwrap_or_default(),
//                     lineage: TaxonomicLineage::default(),
//                     macro_signature: sig.clone(),
//                     meso_signature: sig.clone(),
//                     micro_signature: sig,
//                     levels: vec![],
//                 })
//             })
//             .collect();

//         Ok(ProcessingResult {
//             metrics,
//             signatures: processed_sigs,
//             sample_id,
//         })
//     }

//     pub fn apply_quality_control(
//         &self,
//         sequence: &[u8],
//         min_length: Option<usize>,
//     ) -> Option<Vec<u8>> {
//         let min_len = min_length.unwrap_or(50);
//         if sequence.is_empty() || sequence.len() < min_len {
//             return None;
//         }
//         Some(sequence.to_vec())
//     }

//     pub fn process_chunk(
//         &self,
//         chunk: &[(Vec<u8>, Option<Vec<u8>>)],
//         metrics: &Arc<Mutex<ProcessingMetrics>>,
//         signature: &Arc<Mutex<MultiResolutionSignature>>,
//     ) -> Result<(), ProcessingError> {
//         chunk.par_iter().try_for_each(|(seq, _quality)| {
//             let processed_seq = self.process_sequence(seq)?;
//             if !processed_seq.is_empty() {
//                 // Update metrics
//                 {
//                     let mut metrics = metrics.lock().unwrap();
//                     metrics.processed_reads += 1;
//                     metrics.processed_bases += processed_seq.len();
//                 }

//                 // Update signature at each resolution level
//                 let mut sig_guard = signature.lock().unwrap();
//                 for level in &mut sig_guard.levels {
//                     level.add_sequence(&processed_seq).map_err(|e| {
//                         ProcessingError::SignatureError(format!(
//                             "Signature update failed at level {}: {}",
//                             level.kmer_size, e
//                         ))
//                     })?;
//                 }
//             }
//             Ok(())
//         })
//     }
// }

// impl ProcessingResult {
//     pub fn generate_report(&self) -> String {
//         let mut report = String::new();
//         report.push_str(&format!("Sample ID: {}\n", self.sample_id));
//         report.push_str(&format!("Total reads: {}\n", self.metrics.total_reads));
//         report.push_str(&format!(
//             "Reads passed QC: {}\n",
//             self.metrics.reads_passed_qc
//         ));
//         report.push_str(&format!(
//             "Average read length: {:.2}\n",
//             self.metrics.avg_read_length
//         ));
//         report.push_str(&format!(
//             "Processing time: {:.2}s\n",
//             self.metrics.processing_time_seconds
//         ));
//         report
//     }
// }

// // --- Report Generation --- (Implementation unchanged)
// pub fn generate_report(results: &ClassificationResults) -> Result<String, ProcessingError> {
//     let mut report = String::new();

//     // Add processing metrics
//     report.push_str(&format!("\nProcessing Metrics:\n"));
//     report.push_str(&format!("Total reads: {}\n", results.metrics.total_reads));
//     report.push_str(&format!(
//         "Passed reads: {}\n",
//         results.metrics.reads_passed_qc
//     ));
//     report.push_str(&format!(
//         "Average read length: {:.2}\n",
//         results.metrics.avg_read_length
//     ));
//     report.push_str(&format!(
//         "Processing time: {:.2}s\n",
//         results.metrics.processing_time_seconds
//     ));

//     // Add classification results
//     report.push_str("\nClassification Results:\n");
//     for classification in &results.classifications {
//         report.push_str(&format!(
//             "- Level: {:?}, Taxon: {}\n",
//             classification.level, classification.taxon_id
//         ));
//     }

//     // Add strain abundances if any
//     if !results.strain_abundances.is_empty() {
//         report.push_str("\nStrain Abundances:\n");
//         for (strain_id, (abundance, confidence)) in &results.strain_abundances {
//             report.push_str(&format!(
//                 "- {}: {:.2}% (confidence: {:.2})\n",
//                 strain_id,
//                 abundance * 100.0,
//                 confidence
//             ));
//         }
//     }

//     Ok(report)
// }

// // --- CLI Runner --- (Implementation largely unchanged, emphasizing placeholder nature)
// // TODO: Implement proper CLI argument parsing using `clap`.
// // TODO: Allow configuration of K-mer sizes, sketch size, QC params, etc. via CLI.
// pub fn run_fastq_cli(
//     fastq_path: impl AsRef<Path>,
//     sample_id: &str,
//     db_path: impl AsRef<Path>,
//     output_dir: impl AsRef<Path>,
//     threads: usize,
// ) -> Result<(), ProcessingError> {
//     env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

//     info!("=================================================");
//     info!("===   Starting FASTQ Processing Pipeline    ===");
//     info!("=================================================");
//     // ... (logging inputs) ...

//     // --- Parameter Setup (Hardcoded - Replace with CLI parsing) ---
//     let macro_k = 31;
//     let meso_k = 21;
//     let sketch_size = 1000;
//     let qc_params = QualityControlParams::default();
//     let cache_dir = output_dir.as_ref().join(".cache_sketches");
//     let api_key: Option<String> = None;
//     fs::create_dir_all(&cache_dir).map_err(|e| ProcessingError::IoError(e))?;
//     // ... (logging configuration) ...

//     // --- Initialize Processor ---
//     let init_start = Instant::now();
//     let mut processor = Processor::new();
//     info!(
//         "Processor initialized in {:.2}s.",
//         init_start.elapsed().as_secs_f64()
//     );

//     // --- Run Processing ---
//     info!("Starting file processing for {}...", sample_id);
//     let processing_start = Instant::now();
//     let results = processor.process_file(
//         fastq_path.as_ref().to_path_buf(),
//         sample_id.to_string(),
//         output_dir.as_ref().to_path_buf(),
//     )?;
//     info!(
//         "File processing completed in {:.2}s.",
//         processing_start.elapsed().as_secs_f64()
//     );

//     // --- Generate and Output Report ---
//     info!("Generating final text report...");
//     let report = results.generate_report();
//     println!("\n{}", report); // Print report to stdout

//     // --- Final Confirmation ---
//     if let Some(ref results_json_path) = results.results_file {
//         info!(
//             "Detailed JSON results saved to: {}",
//             results_json_path.display()
//         );
//     } else {
//         warn!("Results file path was not set (e.g., processing might have yielded empty results).");
//     }

//     info!("=================================================");
//     info!("=== Processing Finished Successfully for {} ===", sample_id);
//     info!("=================================================");
//     Ok(())
// }

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use std::path::PathBuf;

//     fn create_test_processor() -> Result<Processor, ProcessingError> {
//         Ok(Processor::new())
//     }

//     #[test]
//     fn test_qc_empty_sequence() {
//         let processor = create_test_processor().unwrap();
//         let empty_seq = vec![];
//         assert!(processor.apply_quality_control(&empty_seq, None).is_none());
//     }

//     #[test]
//     fn test_qc_min_length() {
//         let processor = create_test_processor().unwrap();
//         let short_seq = vec![b'A'; 10]; // Shorter than min length
//         assert!(processor.apply_quality_control(&short_seq, None).is_none());
//     }

//     #[test]
//     fn test_signature_creation() {
//         let processor = create_test_processor().unwrap();
//         let result = processor.create_empty_sample_signature("test_sample");
//         assert!(result.is_ok());
//     }
// }

// // Remove the Hash impl for &[u8] and add a newtype wrapper
// #[derive(Debug)]
// struct ByteSlice<'a>(&'a [u8]);

// impl<'a> Hash for ByteSlice<'a> {
//     fn hash<H: Hasher>(&self, state: &mut H) {
//         self.0.iter().for_each(|byte| byte.hash(state));
//     }
// }
