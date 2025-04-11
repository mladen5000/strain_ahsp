use anyhow::{anyhow, Result};
use bio::io::fastq;
use log::info;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

// Define SequenceRecord with the necessary fields
#[derive(Debug, Clone)]
pub struct SequenceRecord {
    pub id: String,
    pub seq: String,
    pub qual: Option<String>,
}

// Function to find sequence files
fn find_sequence_files(input_paths: &[String]) -> Result<Vec<String>> {
    // Filter for files that exist and have .fastq or .fq extensions
    let mut result = Vec::new();
    for path in input_paths {
        let path_obj = Path::new(path);
        if path_obj.exists() && path_obj.is_file() {
            if let Some(ext) = path_obj.extension() {
                let ext_str = ext.to_string_lossy().to_lowercase();
                if ext_str == "fastq" || ext_str == "fq" {
                    result.push(path.clone());
                }
            }
        }
    }
    Ok(result)
}

pub fn read_sequences_stream(input_paths: &[String]) -> Result<Vec<SequenceRecord>> {
    let files_to_process = find_sequence_files(input_paths)?;
    if files_to_process.is_empty() {
        return Err(anyhow!(
            "No sequence files found in the provided paths: {:?}",
            input_paths
        ));
    }
    info!("Reading sequences from {} file(s).", files_to_process.len());

    // Pre-allocate with a reasonable capacity to avoid frequent reallocations
    let mut all_records = Vec::with_capacity(1000);

    for file_path in files_to_process {
        let file = File::open(&file_path)?;
        let reader = BufReader::with_capacity(1024 * 1024, file); // Use a larger buffer size
        let records = fastq::Reader::new(reader).records();

        for record in records {
            match record {
                Ok(rec) => {
                    let id = rec.id().to_owned();
                    let seq = rec.seq().to_owned();
                    let qual = rec.qual().to_owned();

                    let sequence_record = SequenceRecord {
                        id,
                        seq: String::from_utf8_lossy(&seq).to_string(),
                        qual: Some(String::from_utf8_lossy(&qual).to_string()),
                    };
                    all_records.push(sequence_record);
                }
                Err(e) => {
                    return Err(anyhow!(
                        "Failed to parse record in file {:?}: {}",
                        file_path,
                        e
                    ));
                }
            }
        }
    }

    // Shrink to actual size to minimize memory usage
    all_records.shrink_to_fit();

    info!("Finished reading {} sequences.", all_records.len());
    Ok(all_records)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::tempdir;

    fn create_test_fastq(dir: &Path, name: &str, records: Vec<(String, String, String)>) -> String {
        let file_path = dir.join(name);
        let mut file = File::create(&file_path).unwrap();

        for (id, seq, qual) in records {
            writeln!(file, "@{}", id).unwrap();
            writeln!(file, "{}", seq).unwrap();
            writeln!(file, "+").unwrap();
            writeln!(file, "{}", qual).unwrap();
        }

        file_path.to_str().unwrap().to_string()
    }

    #[test]
    fn test_read_sequences_stream() {
        let temp_dir = tempdir().unwrap();

        // Create a test FASTQ file
        let records = vec![
            ("seq1".to_string(), "ACGT".to_string(), "IIII".to_string()),
            ("seq2".to_string(), "TGCA".to_string(), "HHHH".to_string()),
        ];

        let file_path = create_test_fastq(temp_dir.path(), "test.fastq", records);

        // Test reading the sequences
        let result = read_sequences_stream(&[file_path]).unwrap();

        assert_eq!(result.len(), 2);
        assert_eq!(result[0].id, "seq1");
        assert_eq!(result[0].seq, "ACGT");
        assert_eq!(result[0].qual.as_ref().unwrap(), "IIII");

        assert_eq!(result[1].id, "seq2");
        assert_eq!(result[1].seq, "TGCA");
        assert_eq!(result[1].qual.as_ref().unwrap(), "HHHH");
    }

    #[test]
    fn test_empty_paths() {
        let result = read_sequences_stream(&[]);
        assert!(result.is_err());
    }

    #[test]
    fn test_nonexistent_file() {
        let result = read_sequences_stream(&["nonexistent.fastq".to_string()]);
        assert!(result.is_err());
    }
}
