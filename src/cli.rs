use crate::adaptive::classifier::AdaptiveClassifier;
use crate::database::DatabaseManager;r};
use crate::pipeline::report::Cli as ReportCli; // Changed this line
use crate::visualization::{VisualizationType, Visualizer};
use clap::{Parser, Subcommand};
use std::fs::File;
use std::path::PathBuf;

#[derive(Subcommand)]
pub enum Commands {
    // ... existing commands ...
    /// Generate visualizations for processed results
    Visualize {
        /// Path to the results JSON file
        #[arg(short, long)]
        results: PathBuf,

        /// Path to the output directory
        #[arg(short, long, default_value = "visualizations")]
        output: PathBuf,
    },

    /// Compare multiple samples
    CompareSamples {
        /// Paths to multiple results JSON files
        #[arg(short, long)]
        results: Vec<PathBuf>,

        /// Path to the output directory
        #[arg(short, long, default_value = "comparisons")]
        output: PathBuf,
    },
}

/// Main entry point for CLI
pub fn run_cli(cli: ReportCli) -> Result<(), Box<dyn std::error::Error>> {
    // Updated parameter type
    // ... existing CLI handling ...

    match cli.command {
        // ... existing command handlers ...
        Commands::Visualize { results, output } => {
            // Load results from file
            let file = File::open(&results)?;
            let results_data: ClassificationResults = serde_json::from_reader(file)?;

            // Create visualizer
            let visualizer = Visualizer::new(&output)?;

            // Generate visualizations
            println!(
                "Generating visualizations for sample: {}",
                results_data.sample_id
            );

            // Generate individual visualizations
            let taxonomy_chart = visualizer
                .generate_visualization(&results_data, VisualizationType::TaxonomySunburst)?;
            println!("Generated taxonomy chart: {}", taxonomy_chart.display());

            let strain_chart = visualizer
                .generate_visualization(&results_data, VisualizationType::StrainBarChart)?;
            println!("Generated strain chart: {}", strain_chart.display());

            let confidence_chart = visualizer
                .generate_visualization(&results_data, VisualizationType::ConfidenceHeatmap)?;
            println!("Generated confidence chart: {}", confidence_chart.display());

            // Generate comprehensive HTML report
            let html_report = visualizer.generate_html_report(&results_data)?;
            println!("Generated HTML report: {}", html_report.display());
            println!("Open this file in a web browser to view the interactive report");
        }

        Commands::CompareSamples { results, output } => {
            // Load all result files
            let mut results_data = Vec::new();

            for path in &results {
                let file = File::open(path)?;
                let result: ClassificationResults = serde_json::from_reader(file)?;
                results_data.push(result);
            }

            // Create visualizer
            let visualizer = Visualizer::new(&output)?;

            // Generate comparison visualization
            println!("Comparing {} samples", results_data.len());
            let comparison = visualizer.compare_samples(&results_data)?;
            println!("Generated comparison: {}", comparison.display());
        }
    }

    Ok(())
}
