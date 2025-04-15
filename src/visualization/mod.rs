pub mod cli;
pub mod plotter;

use std::path::{Path, PathBuf};
pub enum VisualizationType {
    TaxonomySunburst,
    StrainBarChart,
    ConfidenceHeatmap,
}

pub struct Visualizer {
    output_dir: PathBuf,
}

impl Visualizer {
    pub fn new(output_dir: &Path) -> Result<Self, std::io::Error> {
        std::fs::create_dir_all(output_dir)?;
        Ok(Self {
            output_dir: output_dir.to_owned(),
        })
    }

    pub fn generate_visualization(
        &self,
        results: &ClassificationResults,
        viz_type: VisualizationType,
    ) -> Result<PathBuf, Box<dyn std::error::Error>> {
        // Implementation needed
        todo!("Implement visualization generation")
    }

    pub fn generate_html_report(
        &self,
        results: &ClassificationResults,
    ) -> Result<PathBuf, Box<dyn std::error::Error>> {
        // Implementation needed
        todo!("Implement HTML report generation")
    }

    pub fn compare_samples(
        &self,
        results: &[ClassificationResults],
    ) -> Result<PathBuf, Box<dyn std::error::Error>> {
        // Implementation needed
        todo!("Implement sample comparison")
    }
}
