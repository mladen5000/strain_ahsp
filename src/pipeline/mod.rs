pub mod processor;
pub mod qc;
pub mod report;

pub use processor::{ClassificationResults, FastqProcessor, ProcessingMetrics};
pub use qc::QualityControlParams;
pub use report::generate_report;
