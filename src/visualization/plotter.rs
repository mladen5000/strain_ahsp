use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{self, Write};
use std::path::{Path, PathBuf};

use plotters::prelude::*;
use serde::Serialize;
use thiserror::Error;

use crate::adaptive::classifier::{Classification, TaxonomicLevel};
use crate::pipeline::processor::ClassificationResults;

#[derive(Error, Debug)]
pub enum VisualizationError {
    #[error("IO error: {0}")]
    IoError(#[from] io::Error),

    #[error("Plot error: {0}")]
    PlotError(String),

    #[error("Template error: {0}")]
    TemplateError(String),
}

/// Visualization types
pub enum VisualizationType {
    /// Taxonomic hierarchy as a sunburst chart
    TaxonomySunburst,

    /// Strain abundances as a bar chart
    StrainBarChart,

    /// Confidence metrics as a heatmap
    ConfidenceHeatmap,

    /// Comparative visualization between samples
    SampleComparison,
}

/// Template for HTML reports
const HTML_TEMPLATE: &str = r#"<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>AHSP Report: {{sample_id}}</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js@3.7.1/dist/chart.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/d3@7"></script>
    <style>
        body {
            font-family: Arial, sans-serif;
            line-height: 1.6;
            color: #333;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
        }
        h1 {
            color: #2c3e50;
            border-bottom: 2px solid #3498db;
            padding-bottom: 10px;
        }
        .metrics {
            background: #f8f9fa;
            padding: 15px;
            border-radius: 5px;
            margin-bottom: 20px;
        }
        .chart-container {
            width: 100%;
            height: 400px;
            margin: 20px 0;
        }
        .flex-container {
            display: flex;
            flex-wrap: wrap;
            justify-content: space-between;
        }
        .flex-item {
            flex: 0 0 48%;
            margin-bottom: 20px;
        }
        table {
            width: 100%;
            border-collapse: collapse;
        }
        th, td {
            padding: 8px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }
        th {
            background-color: #f2f2f2;
        }
        tr:hover {
            background-color: #f5f5f5;
        }
    </style>
</head>
<body>
    <h1>AHSP Analysis Report: {{sample_id}}</h1>
    
    <div class="metrics">
        <h2>Processing Metrics</h2>
        <p><strong>Total reads:</strong> {{total_reads}}</p>
        <p><strong>Passed QC:</strong> {{passed_reads}} ({{qc_percent}}%)</p>
        <p><strong>Average read length:</strong> {{avg_read_length}} bp</p>
        <p><strong>Processing time:</strong> {{processing_time}} seconds</p>
    </div>
    
    <div class="flex-container">
        <div class="flex-item">
            <h2>Taxonomic Classification</h2>
            <div class="chart-container">
                <canvas id="taxonomyChart"></canvas>
            </div>
        </div>
        
        <div class="flex-item">
            <h2>Strain Abundances</h2>
            <div class="chart-container">
                <canvas id="strainChart"></canvas>
            </div>
        </div>
    </div>
    
    <h2>Classification Details</h2>
    <table>
        <tr>
            <th>Rank</th>
            <th>Taxon</th>
            <th>Level</th>
            <th>Confidence</th>
        </tr>
        {{#classifications}}
        <tr>
            <td>{{rank}}</td>
            <td>{{taxon_id}}</td>
            <td>{{level}}</td>
            <td>{{confidence}}</td>
        </tr>
        {{/classifications}}
    </table>
    
    <h2>Strain Details</h2>
    <table>
        <tr>
            <th>Strain</th>
            <th>Abundance</th>
            <th>Confidence Interval</th>
        </tr>
        {{#strains}}
        <tr>
            <td>{{id}}</td>
            <td>{{abundance}}%</td>
            <td>±{{confidence}}%</td>
        </tr>
        {{/strains}}
    </table>
    
    <script>
        // Taxonomy Chart
        const taxonomyCtx = document.getElementById('taxonomyChart').getContext('2d');
        const taxonomyChart = new Chart(taxonomyCtx, {
            type: 'pie',
            data: {
                labels: [{{taxonomy_labels}}],
                datasets: [{
                    data: [{{taxonomy_data}}],
                    backgroundColor: [{{taxonomy_colors}}],
                }]
            },
            options: {
                responsive: true,
                plugins: {
                    legend: {
                        position: 'right',
                    },
                    title: {
                        display: true,
                        text: 'Taxonomic Classification'
                    }
                }
            }
        });
        
        // Strain Chart
        const strainCtx = document.getElementById('strainChart').getContext('2d');
        const strainChart = new Chart(strainCtx, {
            type: 'bar',
            data: {
                labels: [{{strain_labels}}],
                datasets: [{
                    label: 'Abundance (%)',
                    data: [{{strain_data}}],
                    backgroundColor: [{{strain_colors}}],
                    borderColor: [{{strain_border_colors}}],
                    borderWidth: 1,
                    barPercentage: 0.6,
                }]
            },
            options: {
                responsive: true,
                plugins: {
                    legend: {
                        display: false
                    },
                    title: {
                        display: true,
                        text: 'Strain Abundances'
                    }
                },
                scales: {
                    y: {
                        beginAtZero: true,
                        title: {
                            display: true,
                            text: 'Abundance (%)'
                        }
                    }
                }
            }
        });
    </script>
</body>
</html>"#;

/// Visualization generator
pub struct Visualizer {
    /// Output directory for visualizations
    output_dir: PathBuf,
}

impl Visualizer {
    /// Create a new visualizer
    pub fn new(output_dir: impl AsRef<Path>) -> Result<Self, VisualizationError> {
        let output_path = output_dir.as_ref().to_path_buf();

        // Create output directory if it doesn't exist
        if !output_path.exists() {
            fs::create_dir_all(&output_path)?;
        }

        Ok(Visualizer {
            output_dir: output_path,
        })
    }

    /// Generate a visualization from classification results
    pub fn generate_visualization(
        &self,
        results: &ClassificationResults,
        viz_type: VisualizationType,
    ) -> Result<PathBuf, VisualizationError> {
        match viz_type {
            VisualizationType::TaxonomySunburst => self.create_sunburst(results),
            VisualizationType::StrainBarChart => self.create_strain_chart(results),
            VisualizationType::ConfidenceHeatmap => self.create_confidence_heatmap(results),
            VisualizationType::SampleComparison => Err(VisualizationError::PlotError(
                "Sample comparison requires multiple samples".to_string(),
            )),
        }
    }

    /// Generate a comprehensive HTML report
    pub fn generate_html_report(
        &self,
        results: &ClassificationResults,
    ) -> Result<PathBuf, VisualizationError> {
        let output_file = self
            .output_dir
            .join(format!("{}_report.html", results.sample_id));

        // Prepare template data
        let template_data = self.prepare_template_data(results)?;

        // Simple template replacement (in a real implementation, use a templating library)
        let mut html = HTML_TEMPLATE.to_string();

        for (key, value) in template_data {
            html = html.replace(&format!("{{{{{}}}}}", key), &value);
        }

        // Write HTML to file
        let mut file = File::create(&output_file)?;
        file.write_all(html.as_bytes())?;

        Ok(output_file)
    }

    /// Prepare data for HTML template
    fn prepare_template_data(
        &self,
        results: &ClassificationResults,
    ) -> Result<HashMap<String, String>, VisualizationError> {
        let mut data = HashMap::new();

        // Basic information
        data.insert("sample_id".to_string(), results.sample_id.clone());
        data.insert(
            "total_reads".to_string(),
            results.metrics.total_reads.to_string(),
        );
        data.insert(
            "passed_reads".to_string(),
            results.metrics.passed_reads.to_string(),
        );

        let qc_percent = if results.metrics.total_reads > 0 {
            format!(
                "{:.1}",
                100.0 * results.metrics.passed_reads as f64 / results.metrics.total_reads as f64
            )
        } else {
            "0.0".to_string()
        };
        data.insert("qc_percent".to_string(), qc_percent);

        data.insert(
            "avg_read_length".to_string(),
            format!("{:.1}", results.metrics.avg_read_length),
        );
        data.insert(
            "processing_time".to_string(),
            format!("{:.2}", results.metrics.processing_time_seconds),
        );

        // Taxonomic classification data
        let mut taxonomy_labels = String::new();
        let mut taxonomy_data = String::new();
        let mut taxonomy_colors = String::new();

        // Use the first classification as primary
        if let Some(classification) = results.classifications.first() {
            // Extract lineage
            let mut taxa = Vec::new();
            if !classification.lineage.is_empty() {
                for taxon in &classification.lineage {
                    taxa.push(taxon.clone());
                }
            } else {
                taxa.push(classification.taxon_id.clone());
            }

            // Generate chart data
            for (i, taxon) in taxa.iter().enumerate() {
                if i > 0 {
                    taxonomy_labels.push_str(", ");
                    taxonomy_data.push_str(", ");
                    taxonomy_colors.push_str(", ");
                }

                taxonomy_labels.push_str(&format!("'{}'", taxon));
                taxonomy_data.push_str(&format!("{}", 100.0 / taxa.len() as f64)); // Simplified

                // Generate a color based on index
                let hue = (i as f64 * 137.5) % 360.0;
                taxonomy_colors.push_str(&format!("'hsl({}, 70%, 60%)'", hue));
            }
        }

        data.insert("taxonomy_labels".to_string(), taxonomy_labels);
        data.insert("taxonomy_data".to_string(), taxonomy_data);
        data.insert("taxonomy_colors".to_string(), taxonomy_colors);

        // Strain abundance data
        let mut strain_labels = String::new();
        let mut strain_data = String::new();
        let mut strain_colors = String::new();
        let mut strain_border_colors = String::new();

        // Get top strains sorted by abundance
        let mut strains: Vec<_> = results.strain_abundances.iter().collect();
        strains.sort_by(|a, b| b.1 .0.partial_cmp(&a.1 .0).unwrap());

        for (i, (strain_id, (abundance, _))) in strains.iter().enumerate() {
            if i > 0 {
                strain_labels.push_str(", ");
                strain_data.push_str(", ");
                strain_colors.push_str(", ");
                strain_border_colors.push_str(", ");
            }

            strain_labels.push_str(&format!("'{}'", strain_id));
            strain_data.push_str(&format!("{:.2}", abundance * 100.0));

            // Generate a color based on index
            let hue = (i as f64 * 137.5) % 360.0;
            strain_colors.push_str(&format!("'hsla({}, 70%, 60%, 0.7)'", hue));
            strain_border_colors.push_str(&format!("'hsl({}, 70%, 50%)'", hue));
        }

        data.insert("strain_labels".to_string(), strain_labels);
        data.insert("strain_data".to_string(), strain_data);
        data.insert("strain_colors".to_string(), strain_colors);
        data.insert("strain_border_colors".to_string(), strain_border_colors);

        // Classification details
        let mut classifications = String::new();
        for (i, classification) in results.classifications.iter().enumerate() {
            classifications.push_str(&format!(
                "<tr><td>{}</td><td>{}</td><td>{:?}</td><td>{:.2}</td></tr>\n",
                i + 1,
                classification.taxon_id,
                classification.level,
                classification.confidence
            ));
        }
        data.insert("classifications".to_string(), classifications);

        // Strain details
        let mut strain_rows = String::new();
        for (strain_id, (abundance, confidence)) in &results.strain_abundances {
            strain_rows.push_str(&format!(
                "<tr><td>{}</td><td>{:.2}</td><td>±{:.2}</td></tr>\n",
                strain_id,
                abundance * 100.0,
                confidence * 100.0
            ));
        }
        data.insert("strains".to_string(), strain_rows);

        Ok(data)
    }

    /// Create a taxonomy sunburst chart
    fn create_sunburst(
        &self,
        results: &ClassificationResults,
    ) -> Result<PathBuf, VisualizationError> {
        let output_file = self
            .output_dir
            .join(format!("{}_taxonomy_sunburst.svg", results.sample_id));

        // Get dimensions
        let width = 800;
        let height = 600;

        // Create SVG drawing area
        let root = SVGBackend::new(&output_file, (width, height)).into_drawing_area();
        root.fill(&WHITE)?;

        // We'll use a simple pie chart as a substitute for a sunburst
        // (in a real implementation, use D3.js or a more sophisticated library)
        let mut chart = ChartBuilder::on(&root)
            .caption("Taxonomic Classification", ("sans-serif", 30))
            .margin(10)
            .build_cartesian_2d(0.0..1.0, 0.0..1.0)?;

        chart.configure_mesh().disable_mesh().draw()?;

        // Use top classification
        if let Some(classification) = results.classifications.first() {
            let center = (0.5, 0.5);
            let radius = 0.4; // Outer radius

            // Draw pie sections based on lineage
            let mut taxa = Vec::new();
            if !classification.lineage.is_empty() {
                for taxon in &classification.lineage {
                    taxa.push(taxon.clone());
                }
            } else {
                taxa.push(classification.taxon_id.clone());
            }

            let slice_angle = std::f64::consts::PI * 2.0 / taxa.len() as f64;

            for (i, taxon) in taxa.iter().enumerate() {
                let start_angle = i as f64 * slice_angle;
                let end_angle = (i + 1) as f64 * slice_angle;

                // Generate a color based on index
                let hue = (i as f64 * 137.5) % 360.0;
                let color = RGBColor(
                    ((hue + 120.0) % 360.0 / 360.0 * 255.0) as u8,
                    ((hue + 240.0) % 360.0 / 360.0 * 255.0) as u8,
                    (hue / 360.0 * 255.0) as u8,
                );

                // Draw pie slice
                for r in 0..100 {
                    let inner_r = radius * (r as f64 / 100.0);
                    let outer_r = radius * ((r + 1) as f64 / 100.0);

                    root.draw(&Polygon::new(
                        (start_angle..=end_angle)
                            .step(0.1)
                            .map(|angle| {
                                let (sin, cos) = angle.sin_cos();
                                let x = center.0 + outer_r * cos;
                                let y = center.1 + outer_r * sin;
                                (x, y)
                            })
                            .chain((start_angle..=end_angle).step(0.1).rev().map(|angle| {
                                let (sin, cos) = angle.sin_cos();
                                let x = center.0 + inner_r * cos;
                                let y = center.1 + inner_r * sin;
                                (x, y)
                            }))
                            .collect::<Vec<_>>(),
                        &color.mix(0.8 + 0.2 * (r as f64 / 100.0)),
                    ))?;
                }

                // Add label
                let label_angle = start_angle + slice_angle / 2.0;
                let (sin, cos) = label_angle.sin_cos();
                let label_radius = radius * 0.7; // Position label inside the slice
                let label_pos = (center.0 + label_radius * cos, center.1 + label_radius * sin);

                // Add text label
                let style = TextStyle::from(("sans-serif", 12).into_font())
                    .color(&BLACK)
                    .pos(label_pos)
                    .anchor(if cos < 0.0 {
                        TextAlignment::right()
                    } else {
                        TextAlignment::left()
                    });

                let shortened_taxon = if taxon.len() > 15 {
                    format!("{}...", &taxon[0..12])
                } else {
                    taxon.clone()
                };

                root.draw_text(&shortened_taxon, &style)?;
            }
        }

        root.present()?;

        Ok(output_file)
    }

    /// Create a strain abundance bar chart
    fn create_strain_chart(
        &self,
        results: &ClassificationResults,
    ) -> Result<PathBuf, VisualizationError> {
        let output_file = self
            .output_dir
            .join(format!("{}_strain_chart.svg", results.sample_id));

        // Get dimensions
        let width = 800;
        let height = 600;

        // Create SVG drawing area
        let root = SVGBackend::new(&output_file, (width, height)).into_drawing_area();
        root.fill(&WHITE)?;

        // Get top strains sorted by abundance
        let mut strains: Vec<_> = results.strain_abundances.iter().collect();
        strains.sort_by(|a, b| b.1 .0.partial_cmp(&a.1 .0).unwrap());

        // Limit to top 10 strains for clarity
        let top_strains: Vec<_> = strains.into_iter().take(10).collect();

        if top_strains.is_empty() {
            // No strain data, draw empty chart with message
            let mut chart = ChartBuilder::on(&root)
                .caption("Strain Abundances", ("sans-serif", 30))
                .margin(10)
                .build_cartesian_2d(0..1, 0..1)?;

            chart.configure_mesh().disable_mesh().draw()?;

            // Add message
            root.draw_text(
                "No strain abundance data available",
                &TextStyle::from(("sans-serif", 20).into_font())
                    .color(&BLACK)
                    .pos(0.5, 0.5)
                    .anchor(TextAlignment::center()),
            )?;
        } else {
            // Create strain labels and abundances
            let labels: Vec<String> = top_strains
                .iter()
                .map(|(id, _)| {
                    if id.len() > 15 {
                        format!("{}...", &id[0..12])
                    } else {
                        id.clone()
                    }
                })
                .collect();

            let values: Vec<f64> = top_strains
                .iter()
                .map(|(_, (abundance, _))| abundance * 100.0)
                .collect();

            // Find max value for scaling
            let max_value = values.iter().cloned().fold(0.0, f64::max).max(5.0);

            // Create bar chart
            let mut chart = ChartBuilder::on(&root)
                .caption("Strain Abundances", ("sans-serif", 30))
                .margin(50)
                .x_label_area_size(40)
                .y_label_area_size(60)
                .build_cartesian_2d(0..labels.len(), 0.0..max_value)?;

            chart
                .configure_mesh()
                .disable_x_mesh()
                .y_desc("Abundance (%)")
                .draw()?;

            // Draw bars
            let mut bar_series = Vec::new();

            for (i, &value) in values.iter().enumerate() {
                // Generate a color based on index
                let hue = (i as f64 * 137.5) % 360.0;
                let color = RGBColor(
                    ((hue + 120.0) % 360.0 / 360.0 * 255.0) as u8,
                    ((hue + 240.0) % 360.0 / 360.0 * 255.0) as u8,
                    (hue / 360.0 * 255.0) as u8,
                );

                bar_series.push((i, value, color));
            }

            // Draw the bars
            chart.draw_series(bar_series.iter().map(|&(i, value, color)| {
                let mut bar = Rectangle::new([(i, 0.0), (i + 1, value)], color.mix(0.8).filled());
                bar.set_margin(0, 0, 5, 5);
                bar
            }))?;

            // Add labels
            for (i, label) in labels.iter().enumerate() {
                chart.draw_series(std::iter::once(Text::new(
                    label.clone(),
                    (i + 0.5, -0.5),
                    ("sans-serif", 15.0).into_font().color(&BLACK).rotate(45.0),
                )))?;
            }

            // Add values on top of bars
            for (i, &value) in values.iter().enumerate() {
                chart.draw_series(std::iter::once(Text::new(
                    format!("{:.1}%", value),
                    (i + 0.5, value + 0.5),
                    ("sans-serif", 15.0).into_font().color(&BLACK),
                )))?;
            }
        }

        root.present()?;

        Ok(output_file)
    }

    /// Create a confidence heatmap
    fn create_confidence_heatmap(
        &self,
        results: &ClassificationResults,
    ) -> Result<PathBuf, VisualizationError> {
        let output_file = self
            .output_dir
            .join(format!("{}_confidence_heatmap.svg", results.sample_id));

        // Get dimensions
        let width = 800;
        let height = 600;

        // Create SVG drawing area
        let root = SVGBackend::new(&output_file, (width, height)).into_drawing_area();
        root.fill(&WHITE)?;

        // Collect confidence values at different taxonomic levels
        let mut confidence_data = Vec::new();

        for classification in &results.classifications {
            let level = match classification.level {
                TaxonomicLevel::Domain => "Domain",
                TaxonomicLevel::Phylum => "Phylum",
                TaxonomicLevel::Class => "Class",
                TaxonomicLevel::Order => "Order",
                TaxonomicLevel::Family => "Family",
                TaxonomicLevel::Genus => "Genus",
                TaxonomicLevel::Species => "Species",
                TaxonomicLevel::StrainGroup => "Strain Group",
                TaxonomicLevel::Strain => "Strain",
                TaxonomicLevel::Unknown => "Unknown",
            };

            confidence_data.push((level.to_string(), classification.confidence));
        }

        // Add strain confidence if available
        for (strain_id, (_, confidence)) in results.strain_abundances.iter().take(5) {
            confidence_data.push((
                format!(
                    "Strain: {}",
                    if strain_id.len() > 12 {
                        format!("{}...", &strain_id[0..9])
                    } else {
                        strain_id.clone()
                    }
                ),
                *confidence,
            ));
        }

        // Sort by confidence descending
        confidence_data.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

        if confidence_data.is_empty() {
            // No confidence data available
            let mut chart = ChartBuilder::on(&root)
                .caption("Confidence Metrics", ("sans-serif", 30))
                .margin(10)
                .build_cartesian_2d(0..1, 0..1)?;

            chart.configure_mesh().disable_mesh().draw()?;

            // Add message
            root.draw_text(
                "No confidence data available",
                &TextStyle::from(("sans-serif", 20).into_font())
                    .color(&BLACK)
                    .pos(0.5, 0.5)
                    .anchor(TextAlignment::center()),
            )?;
        } else {
            // Create labels and values
            let labels: Vec<String> = confidence_data
                .iter()
                .map(|(level, _)| level.clone())
                .collect();
            let values: Vec<f64> = confidence_data.iter().map(|(_, conf)| *conf).collect();

            // Create horizontal bar chart
            let mut chart = ChartBuilder::on(&root)
                .caption("Confidence Metrics", ("sans-serif", 30))
                .margin(5)
                .x_label_area_size(40)
                .y_label_area_size(120)
                .build_cartesian_2d(0.0..1.0, 0..labels.len())?;

            chart
                .configure_mesh()
                .disable_y_mesh()
                .x_desc("Confidence")
                .y_desc("Taxonomic Level")
                .draw()?;

            // Draw bars
            chart.draw_series(values.iter().enumerate().map(|(i, &value)| {
                // Color based on confidence
                let color = if value < 0.5 {
                    RGBColor(255, 50, 50)
                } else if value < 0.8 {
                    RGBColor(255, 200, 50)
                } else {
                    RGBColor(50, 200, 50)
                };

                let mut bar = Rectangle::new([(0.0, i), (value, i + 1)], color.mix(0.8).filled());
                bar.set_margin(5, 5, 5, 5);
                bar
            }))?;

            // Add labels
            for (i, label) in labels.iter().enumerate() {
                chart.draw_series(std::iter::once(Text::new(
                    label.clone(),
                    (-0.05, i as f64 + 0.5),
                    ("sans-serif", 15.0).into_font().color(&BLACK),
                )))?;
            }

            // Add values at end of bars
            for (i, &value) in values.iter().enumerate() {
                chart.draw_series(std::iter::once(Text::new(
                    format!("{:.2}", value),
                    (value + 0.02, i as f64 + 0.5),
                    ("sans-serif", 15.0).into_font().color(&BLACK),
                )))?;
            }
        }

        root.present()?;

        Ok(output_file)
    }

    /// Compare multiple samples (stub for future implementation)
    pub fn compare_samples(
        &self,
        results: &[ClassificationResults],
    ) -> Result<PathBuf, VisualizationError> {
        if results.len() < 2 {
            return Err(VisualizationError::PlotError(
                "At least two samples are required for comparison".to_string(),
            ));
        }

        // This would be implemented with more sophisticated visualization
        // For now, we'll create a simple bar chart comparing the top strains

        let output_file = self.output_dir.join("sample_comparison.svg");

        // Get dimensions
        let width = 1000;
        let height = 600;

        // Create SVG drawing area
        let root = SVGBackend::new(&output_file, (width, height)).into_drawing_area();
        root.fill(&WHITE)?;

        // Draw a placeholder message for now
        let mut chart = ChartBuilder::on(&root)
            .caption("Sample Comparison", ("sans-serif", 30))
            .margin(10)
            .build_cartesian_2d(0..1, 0..1)?;

        chart.configure_mesh().disable_mesh().draw()?;

        // Add message
        root.draw_text(
            &format!(
                "Comparison of {} samples (to be implemented)",
                results.len()
            ),
            &TextStyle::from(("sans-serif", 20).into_font())
                .color(&BLACK)
                .pos(0.5, 0.5)
                .anchor(TextAlignment::center()),
        )?;

        root.present()?;

        Ok(output_file)
    }
}

/// Add visualization capability to FASTQ processor
impl crate::pipeline::processor::FastqProcessor {
    /// Generate visualizations for a processed sample
    pub fn generate_visualizations(
        &self,
        results: &ClassificationResults,
        output_dir: impl AsRef<Path>,
    ) -> Result<Vec<PathBuf>, VisualizationError> {
        let visualizer = Visualizer::new(output_dir)?;

        // Generate individual visualizations
        let mut output_files = Vec::new();

        // Taxonomy sunburst
        if let Ok(path) =
            visualizer.generate_visualization(results, VisualizationType::TaxonomySunburst)
        {
            output_files.push(path);
        }

        // Strain bar chart
        if let Ok(path) =
            visualizer.generate_visualization(results, VisualizationType::StrainBarChart)
        {
            output_files.push(path);
        }

        // Confidence heatmap
        if let Ok(path) =
            visualizer.generate_visualization(results, VisualizationType::ConfidenceHeatmap)
        {
            output_files.push(path);
        }

        // HTML report (comprehensive)
        if let Ok(path) = visualizer.generate_html_report(results) {
            output_files.push(path);
        }

        Ok(output_files)
    }
}
