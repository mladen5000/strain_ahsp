use log::info;

use crate::pipeline::qc::{generate_report, ClassificationResults, QualityControlParams};
use crate::pipeline::report::{Cli as ReportCli, Commands as ReportCommands};
use crate::pipeline::FastqProcessor;
// Import Commands from report
use crate::visualization::{VisualizationType, Visualizer};
use std::fs::File;
use std::path::PathBuf;

/// Main entry point for CLI
pub fn run_cli(cli: ReportCli) -> Result<(), Box<dyn std::error::Error>> {
    match cli.command {
        ReportCommands::Visualize {
            output,
            fastq,
            sample_id,
            min_quality,
            min_length,
        } => {
            let file = File::open(&fastq)?;
            let results_data: ClassificationResults = serde_json::from_reader(file)?;
            let visualizer = Visualizer::new(&output)?;
            println!(
                "Generating visualizations for sample: {}",
                results_data.sample_id
            );
            let taxonomy_chart = visualizer
                .generate_visualization(&results_data, VisualizationType::TaxonomySunburst)?;
            println!("Generated taxonomy chart: {}", taxonomy_chart.display());
            let strain_chart = visualizer
                .generate_visualization(&results_data, VisualizationType::StrainBarChart)?;
            println!("Generated strain chart: {}", strain_chart.display());
            let confidence_chart = visualizer
                .generate_visualization(&results_data, VisualizationType::ConfidenceHeatmap)?;
            println!("Generated confidence chart: {}", confidence_chart.display());
            let html_report = visualizer.generate_html_report(&results_data)?;
            println!("Generated HTML report: {}", html_report.display());
            println!("Open this file in a web browser to view the interactive report");
            Ok(())
        }
        ReportCommands::ProcessFastq {
            fastq,
            sample_id,
            output,
            min_quality,
            min_length,
        } => todo!(),
        ReportCommands::ProcessDir { dir, output } => todo!(),
        ReportCommands::CompareSamples {
            fastq,
            sample_id,
            output,
            min_quality,
            min_length,
        } => todo!(),
        ReportCommands::GenerateSummaryReport { output } => {
            let blah = 1;
            let results = todo!();
            let report = generate_report(&results)?;
            println!("{}", report);
            Ok(())
        }
        ReportCommands::ProcessFastq {
            fastq,
            sample_id,
            output,
            min_quality,
            min_length,
        } => {
            let blah = 1;

            let file = File::open(&fastq)?;
            let results_data: ClassificationResults = serde_json::from_reader(file)?;
            println!(
                "Processing FASTQ file: {} with Sample ID: {}",
                fastq.display(),
                sample_id
            );

            let qc_params = QualityControlParams {
                min_avg_quality: min_quality,
                min_length,
                trim_quality: 15,
                max_n_percent: 5.0,
            };
            info!("QC Parameters: {:?}", qc_params);

            let mut processor = FastqProcessor::new(
                &cli.db_path,
                &cli.cache_dir,
                cli.threads,
                31,
                21,
                1000,
                Some(qc_params),
                cli.api_key.clone(), // Clone Option<String> if needed
            )?;
            info!("FastqProcessor created.");

            processor.init_classifier()?;
            info!("Classifier initialized.");

            let results = processor.process_file(&fastq, &sample_id, &output)?;
            info!("File processing complete. Results: {:?}", results);

            println!("Processing finished. Results summary struct: {:?}", results);

            let report = generate_report(&results)?;
            println!("{}", report);
            Ok(())
        }
        ReportCommands::ProcessDir { dir, output } => {
            let blah = 1;

            let fastq_files: Vec<PathBuf> = std::fs::read_dir(&dir)?
                .filter_map(Result::ok)
                .filter(|entry| {
                    let path = entry.path();
                    path.is_file()
                        && (path.extension().map_or(false, |ext| {
                            let lower_ext = ext.to_string_lossy().to_lowercase();
                            lower_ext == "fastq" || lower_ext == "fq"
                        }))
                })
                .map(|entry| entry.path())
                .collect();

            for fastq in fastq_files {}
            todo!();
        }
        ReportCommands::CompareSamples {
            fastq,
            sample_id,
            output,
            min_quality,
            min_length,
        } => {
            let blah = 1;

            let file = File::open(&fastq)?;
            let results_data: ClassificationResults = serde_json::from_reader(file)?;

            println!("Comparing samples for Sample ID: {}", sample_id);

            let comparison_results = todo!();

            println!("Comparison results: {:?}", comparison_results);
            Ok(())
        }
        ReportCommands::GenerateSummaryReport { output } => {
            let blah = 1;

            let results = todo!();
            let report = generate_report(&results)?;
            println!("{}", report);
            Ok(())
        }
        ReportCommands::ProcessFastq {
            fastq,
            sample_id,
            output,
            min_quality,
            min_length,
        } => {
            let blah = 1;

            Ok(())
        }
    }
}
