use log::{debug, error};
use rayon::prelude::*;
use std::sync::Arc;
use std::thread;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum ParallelError {
    #[error("Thread error: {0}")]
    ThreadError(String),

    #[error("Processing error: {0}")]
    ProcessingError(String),

    #[error("Invalid chunk size: {0}")]
    InvalidChunkSize(usize),
}

/// Configuration for parallel processing
pub struct ParallelConfig {
    /// Number of threads to use
    pub threads: usize,

    /// Size of chunks for processing
    pub chunk_size: usize,

    /// Whether to continue on errors
    pub continue_on_error: bool,
}

impl Default for ParallelConfig {
    fn default() -> Self {
        ParallelConfig {
            threads: num_cpus::get(),
            chunk_size: 1000,
            continue_on_error: false,
        }
    }
}

/// Process items in parallel using rayon
pub fn parallel_process<T, U, F, E>(
    items: Vec<T>,
    processor: F,
    config: Option<ParallelConfig>,
) -> Result<Vec<U>, E>
where
    T: Send + Sync,
    U: Send,
    F: Fn(&T) -> Result<U, E> + Send + Sync,
    E: From<ParallelError> + Send + Sync,
{
    let config = config.unwrap_or_default();

    // Validate configuration
    if config.chunk_size == 0 {
        return Err(ParallelError::InvalidChunkSize(0).into());
    }

    // Use rayon's thread pool with specified number of threads
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(config.threads)
        .build()
        .map_err(|e| ParallelError::ThreadError(format!("Failed to build thread pool: {}", e)))?;

    // Process in parallel
    pool.install(|| {
        let results: Vec<Result<U, E>> = items
            .par_chunks(config.chunk_size)
            .flat_map(|chunk| {
                chunk
                    .par_iter()
                    .map(|item| {
                        let result = processor(item);
                        if result.is_err() && !config.continue_on_error {
                            error!("Processing error, stopping due to continue_on_error=false");
                        }
                        result
                    })
                    .collect::<Vec<_>>()
            })
            .collect();

        // Handle errors
        if !config.continue_on_error && results.iter().any(|r| r.is_err()) {
            // Return the first error
            for result in results {
                if let Err(e) = result {
                    return Err(e);
                }
            }
            // This should not be reached, but just in case
            Err(ParallelError::ProcessingError(
                "Unknown error during parallel processing".to_string(),
            )
            .into())
        } else {
            // Filter out errors if continue_on_error is true
            let successful: Vec<U> = results
                .into_iter()
                .filter_map(|r| match r {
                    Ok(result) => Some(result),
                    Err(e) => {
                        error!("Error during parallel processing: {:?}", e);
                        None
                    }
                })
                .collect();

            Ok(successful)
        }
    })
}

/// Process items in batches with automatic work distribution
pub fn process_in_batches<T, U, F, E>(
    items: Vec<T>,
    processor: F,
    batch_size: usize,
) -> Result<Vec<U>, E>
where
    T: Send + Sync,
    U: Send,
    F: Fn(Vec<T>) -> Result<Vec<U>, E> + Send + Sync + Clone,
    E: From<ParallelError> + Send + Sync,
{
    if batch_size == 0 {
        return Err(ParallelError::InvalidChunkSize(0).into());
    }

    let num_batches = (items.len() + batch_size - 1) / batch_size;
    debug!(
        "Processing {} items in {} batches of size {}",
        items.len(),
        num_batches,
        batch_size
    );

    if num_batches == 0 {
        return Ok(Vec::new());
    }

    // Create batches
    let batches: Vec<Vec<T>> = (0..num_batches)
        .map(|i| {
            let start = i * batch_size;
            let end = std::cmp::min((i + 1) * batch_size, items.len());
            items[start..end].to_vec()
        })
        .collect();

    // Process batches in parallel
    let processor_arc = Arc::new(processor);
    let mut handles = Vec::with_capacity(num_batches);

    for batch in batches {
        let processor_clone = Arc::clone(&processor_arc);
        let handle = thread::spawn(move || processor_clone(batch));
        handles.push(handle);
    }

    // Collect results
    let mut results = Vec::new();
    for handle in handles {
        match handle.join() {
            Ok(batch_result) => match batch_result {
                Ok(batch_items) => results.extend(batch_items),
                Err(e) => return Err(e),
            },
            Err(_) => {
                return Err(ParallelError::ThreadError("Thread panicked".to_string()).into());
            }
        }
    }

    Ok(results)
}

/// Execute a function with progress tracking on a thread pool
pub struct ParallelExecutor {
    /// Thread pool
    pool: rayon::ThreadPool,

    /// Configuration
    config: ParallelConfig,
}

impl ParallelExecutor {
    /// Create a new parallel executor
    pub fn new(config: Option<ParallelConfig>) -> Result<Self, ParallelError> {
        let config = config.unwrap_or_default();

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(config.threads)
            .build()
            .map_err(|e| {
                ParallelError::ThreadError(format!("Failed to build thread pool: {}", e))
            })?;

        Ok(ParallelExecutor { pool, config })
    }

    /// Execute a function in parallel
    pub fn execute<T, U, F, E>(&self, items: Vec<T>, processor: F) -> Result<Vec<U>, E>
    where
        T: Send + Sync,
        U: Send,
        F: Fn(&T) -> Result<U, E> + Send + Sync,
        E: From<ParallelError> + Send + Sync,
    {
        self.pool
            .install(|| parallel_process(items, processor, Some(self.config.clone())))
    }
}
