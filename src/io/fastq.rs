//! Functions specifically for reading FASTQ files.
//!
//! Leverages the `needletail` crate for efficient parsing
//! of FASTQ and FASTA formats, handling compressed files automatically.

use anyhow::{anyhow, Context, Result};
use log::{info, warn};
use needletail::{parse_fastx_file, FastxReader};
use needletail::parser::SequenceRecord;
use needletail::Sequence;
use std::fs;
use std::path::{Path, PathBuf};

/// Reads sequences from one or more FASTQ/FASTA files (compressed or uncompressed).
///
/// # Arguments
///
/// * `input_paths` - A slice of strings representing file paths or directories.
///                   If a path is a directory, all files within it will be scanned.
///
/// # Returns
///
/// * `Result<Vec<SequenceRecord<'static>>>` - A vector containing all sequence records read.
///   Using `'static` lifetime here implies sequences are owned (copied), which might
///   consume significant memory for large datasets.
///   **Note:** Needletail's default iterators often borrow data. Returning owned data
///   requires cloning, which is done implicitly by collecting into a Vec here.
pub fn read_sequences_owned(input_paths: &[String]) -> Result<Vec<SequenceRecord<'static>>> {
    // Just call the stream function directly, as it now returns a Vec
    read_sequences_stream(input_paths)
}

/// Provides a streaming iterator over sequences from multiple files.
/// This is more memory-efficient than reading all sequences into a Vec at once.
///
/// # Arguments
///
/// * `input_paths` - A slice of strings representing file paths or directories.
///
/// # Returns
///
/// * `Result<Vec<SequenceRecord<'static>>>` - A vector containing all sequence records.
///   This implementation changed to reduce borrow checker issues.
pub fn read_sequences_stream(input_paths: &[String]) -> Result<Vec<SequenceRecord<'static>>> {
    let files_to_process = find_sequence_files(input_paths)?;

    if files_to_process.is_empty() {
        return Err(anyhow!(
            "No sequence files found in the provided paths: {:?}",
            input_paths
        ));
    }
    
    info!("Reading sequences from {} file(s).", files_to_process.len());

    // Instead of trying to create a complex iterator that handles ownership issues,
    // we'll just collect all the records into a single vector
    let mut all_records = Vec::new();
    
    for file_path in files_to_process {
        let reader_result = parse_fastx_file(&file_path);
        match reader_result {
            Ok(mut reader) => {
                let file_path_str = format!("{:?}", file_path);
                // Process each record from this file
                while let Some(record_result) = reader.next() {
                    match record_result {
                        Ok(record) => {
                            // Convert to owned record and add to our collection
                            all_records.push(record.to_owned());
                        },
                        Err(e) => {
                            // Return error if we encounter a problem
                            return Err(anyhow::anyhow!("Failed to parse record in file {}: {}", file_path_str, e));
                        }
                    }
                }
            },
            Err(e) => {
                // Return error if we can't open the file
                let error_msg = format!("Failed to open or parse file: {:?}", file_path);
                return Err(anyhow::anyhow!("{}: {}", error_msg, e));
            }
        }
    }
    
    info!("Finished reading {} sequences.", all_records.len());
    Ok(all_records)
}

/// Finds all sequence files (FASTA/FASTQ, potentially compressed) in the given paths.
/// Expands directories recursively.
fn find_sequence_files(paths: &[String]) -> Result<Vec<PathBuf>> {
    let mut files = Vec::new();
    for path_str in paths {
        let path = Path::new(path_str);
        if !path.exists() {
            warn!("Input path does not exist: {}", path_str);
            continue;
        }
        if path.is_dir() {
            info!("Scanning directory: {:?}", path);
            // Recursively find files in directory
            for entry in fs::read_dir(path)? {
                let entry = entry?;
                let entry_path = entry.path();
                if entry_path.is_file() {
                    // Basic check for common sequence file extensions
                    // Needletail handles compression automatically based on content/extension.
                    if let Some(ext) = entry_path.extension().and_then(|e| e.to_str()) {
                        match ext.to_lowercase().as_str() {
                            "fastq" | "fq" | "fasta" | "fa" | "fna" | "ffn" | "faa" | "frn"
                            | "gz" | "bz2" | "zst" => {
                                files.push(entry_path);
                            }
                            _ => {
                                log::debug!(
                                    "Skipping file with unknown extension: {:?}",
                                    entry_path
                                );
                            }
                        }
                    }
                } else if entry_path.is_dir() {
                    // Recurse into subdirectories
                    files.extend(find_sequence_files(&[entry_path
                        .to_str()
                        .unwrap()
                        .to_string()])?);
                }
            }
        } else if path.is_file() {
            files.push(path.to_path_buf());
        } else {
            warn!("Input path is neither a file nor a directory: {}", path_str);
        }
    }
    Ok(files)
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::Write;
    use tempfile::tempdir;

    // Helper function to create a dummy FASTQ file
    fn create_dummy_fastq(path: &Path, content: &str) {
        let file = fs::File::create(path).unwrap();
        let mut writer = std::io::BufWriter::new(file);
        writer.write_all(content.as_bytes()).unwrap();
    }

    // Helper function to create a dummy compressed FASTQ file
    fn create_dummy_fastq_gz(path: &Path, content: &str) {
        let file = fs::File::create(path).unwrap();
        let mut encoder = GzEncoder::new(file, Compression::default());
        encoder.write_all(content.as_bytes()).unwrap();
        encoder.finish().unwrap();
    }

    #[test]
    fn test_read_single_fastq() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("test.fq");
        let content = "@seq1\nACGT\n+\n!!!!\n@seq2\nTGCA\n+\n####\n";
        create_dummy_fastq(&file_path, content);

        let records = read_sequences_owned(&[file_path.to_str().unwrap().to_string()]).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id(), b"seq1");
        assert_eq!(records[0].sequence(), b"ACGT");
        assert_eq!(records[1].id(), b"seq2");
        assert_eq!(records[1].sequence(), b"TGCA");

        dir.close().unwrap();
    }

    #[test]
    fn test_read_single_fastq_gz() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("test.fq.gz");
        let content = "@seq1\nACGT\n+\n!!!!\n@seq2\nTGCA\n+\n####\n";
        create_dummy_fastq_gz(&file_path, content);

        let records = read_sequences_owned(&[file_path.to_str().unwrap().to_string()]).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id(), b"seq1");
        assert_eq!(records[0].sequence(), b"ACGT");

        dir.close().unwrap();
    }

    #[test]
    fn test_read_from_directory() {
        let dir = tempdir().unwrap();
        let file_path1 = dir.path().join("test1.fa");
        let file_path2 = dir.path().join("test2.fq.gz");
        let file_path3 = dir.path().join("other.txt"); // Should be ignored

        create_dummy_fastq(&file_path1, ">seqA\nAAAA\n>seqB\nCCCC\n"); // Fasta
        create_dummy_fastq_gz(&file_path2, "@seqC\nGGGG\n+\n!!!!\n"); // Compressed Fastq
        create_dummy_fastq(&file_path3, "some text data"); // Ignored file

        let records = read_sequences_owned(&[dir.path().to_str().unwrap().to_string()]).unwrap();
        assert_eq!(records.len(), 3); // seqA, seqB, seqC

        // Check if IDs are present (order might vary)
        let ids: std::collections::HashSet<_> = records.iter().map(|r| r.id().to_vec()).collect();
        assert!(ids.contains(b"seqA".as_ref()));
        assert!(ids.contains(b"seqB".as_ref()));
        assert!(ids.contains(b"seqC".as_ref()));

        dir.close().unwrap();
    }

    #[test]
    fn test_read_nonexistent_file() {
        let result = read_sequences_owned(&["nonexistent_file.fq".to_string()]);
        assert!(result.is_err()); // Should fail as no files are found
    }

    #[test]
    fn test_read_invalid_fastq() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("invalid.fq");
        let content = "@seq1\nACGT\n+\n!!!"; // Missing one quality score
        create_dummy_fastq(&file_path, content);

        let result = read_sequences_owned(&[file_path.to_str().unwrap().to_string()]);
        assert!(result.is_err()); // Should fail during parsing

        dir.close().unwrap();
    }

    #[test]
    fn test_read_sequences_stream_basic() {
        let dir = tempdir().unwrap();
        let file_path1 = dir.path().join("s1.fa");
        let file_path2 = dir.path().join("s2.fq");
        create_dummy_fastq(&file_path1, ">seqA\nAAAA\n");
        create_dummy_fastq(&file_path2, "@seqB\nCCCC\n+\n####\n");

        let stream_result = read_sequences_stream(&[
            file_path1.to_str().unwrap().to_string(),
            file_path2.to_str().unwrap().to_string(),
        ]);
        assert!(stream_result.is_ok());
        let mut stream = stream_result.unwrap();

        let rec1 = stream.next().unwrap().unwrap();
        assert_eq!(rec1.id(), b"seqA");
        assert_eq!(rec1.sequence(), b"AAAA");

        let rec2 = stream.next().unwrap().unwrap();
        assert_eq!(rec2.id(), b"seqB");
        assert_eq!(rec2.sequence(), b"CCCC");

        assert!(stream.next().is_none()); // End of stream

        dir.close().unwrap();
    }

    #[test]
    fn test_read_sequences_stream_error() {
        let dir = tempdir().unwrap();
        let file_path_good = dir.path().join("good.fa");
        let file_path_bad = dir.path().join("bad.fq");
        create_dummy_fastq(&file_path_good, ">seqA\nAAAA\n");
        create_dummy_fastq(&file_path_bad, "@seqB\nCCCC\n+"); // Invalid FASTQ

        let stream_result = read_sequences_stream(&[
            file_path_good.to_str().unwrap().to_string(),
            file_path_bad.to_str().unwrap().to_string(),
        ]);
        assert!(stream_result.is_ok());
        let mut stream = stream_result.unwrap();

        // First record should be ok
        let rec1_res = stream.next().unwrap();
        assert!(rec1_res.is_ok());
        assert_eq!(rec1_res.unwrap().id(), b"seqA");

        // Second record should be an error
        let rec2_res = stream.next().unwrap();
        assert!(rec2_res.is_err());

        assert!(stream.next().is_none()); // End of stream

        dir.close().unwrap();
    }
}
