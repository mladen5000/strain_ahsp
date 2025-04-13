use log::{debug, error};
use rayon::prelude::*;
use std::fmt::Debug; // Import Debug for logging errors
use std::sync::Arc;
use std::thread;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum ParallelError {
    #[error("Thread error: {0}")]
    ThreadError(String),

    #[error("Processing error: {0}")]
    ProcessingError(String),

    #[error("Invalid chunk/batch size: {0}. Must be greater than 0.")]
    InvalidChunkSize(usize),

    #[error("Failed to build thread pool: {0}")]
    ThreadPoolBuildError(String),
}

/// Configuration for parallel processing
#[derive(Clone, Debug)]
pub struct ParallelConfig {
    /// Number of threads to use
    pub threads: usize,

    /// Size of chunks for rayon processing (if applicable, less critical with par_iter)
    pub chunk_size: usize,

    /// Whether to continue processing remaining items if one item fails
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

/// Process items in parallel using rayon's default global pool or a specified pool.
///
/// This function processes items individually in parallel.
/// It handles errors based on the `continue_on_error` flag in the config.
/// If `continue_on_error` is false (default), it returns the first error encountered.
/// If `continue_on_error` is true, it processes all items, logs errors, and returns a Vec
/// containing only the successful results.
pub fn parallel_process<T, U, F, E>(
    items: Vec<T>,
    processor: F,
    config: Option<ParallelConfig>,
) -> Result<Vec<U>, E>
where
    T: Send + Sync,
    U: Send,
    F: Fn(&T) -> Result<U, E> + Send + Sync,
    E: From<ParallelError> + Send + Debug,
{
    let config = config.unwrap_or_default();

    let results: Vec<Result<U, E>> = items.par_iter().map(|item| processor(item)).collect();

    let mut successful_results = Vec::new();
    let mut errors: Vec<E> = Vec::new();

    for result in results {
        match result {
            Ok(value) => {
                successful_results.push(value); // Always collect successes initially
            }
            Err(e) => {
                errors.push(e);
            }
        }
    }

    if !errors.is_empty() {
        if config.continue_on_error {
            for e in errors {
                error!("Error during parallel processing (ignored): {:?}", e);
            }
            // Filter successful results only if we continued on error and errors occurred
            // This seems counter-intuitive; successful_results already *only* contains successes.
            // The filtering logic should be applied *before* deciding whether to return Ok or Err.
            // Let's return the successes collected so far.
            Ok(successful_results)
        } else {
            error!(
                "Error during parallel processing (stopping): {:?}",
                errors[0]
            );
            Err(errors.remove(0)) // Return first error
        }
    } else {
        // No errors, return all collected results (which must be successes)
        Ok(successful_results)
    }
}

/// Process items in batches using standard library threads.
///
/// Each batch is processed by a separate thread. The processor function `F`
/// must be able to handle an owned `Vec<T>` for a batch.
///
/// **Fix:** Added `'static` lifetime bounds to `T`, `U`, `F`, and `E` because
/// `std::thread::spawn` requires the closure and all captured data to live
/// for the static lifetime.
pub fn process_in_batches<T, U, F, E>(
    items: Vec<T>,
    processor: F,
    batch_size: usize,
) -> Result<Vec<U>, E>
where
    T: Send + Sync + Clone + 'static, // Added 'static
    U: Send + 'static,                // Added 'static
    F: Fn(Vec<T>) -> Result<Vec<U>, E> + Send + Sync + Clone + 'static, // Added 'static
    E: From<ParallelError> + Send + Debug + 'static, // Added 'static
{
    if batch_size == 0 {
        return Err(ParallelError::InvalidChunkSize(0).into());
    }

    let item_count = items.len();
    if item_count == 0 {
        return Ok(Vec::new());
    }

    let num_batches = (item_count + batch_size - 1) / batch_size;
    debug!(
        "Processing {} items in {} batches of size {}",
        item_count, num_batches, batch_size
    );

    let batches: Vec<Vec<T>> = items
        .chunks(batch_size)
        .map(|slice| slice.to_vec())
        .collect();

    let processor_arc = Arc::new(processor);
    let mut handles = Vec::with_capacity(num_batches);

    for batch in batches {
        let processor_clone = Arc::clone(&processor_arc);
        // `thread::spawn` requires the closure and its captures ('processor_clone', 'batch')
        // to be 'static. This is why the 'static bounds were added to T, U, F, E.
        let handle = thread::spawn(move || processor_clone(batch));
        handles.push(handle);
    }

    let mut results = Vec::with_capacity(item_count); // Pre-allocate roughly
    for handle in handles {
        match handle.join() {
            Ok(batch_result) => {
                match batch_result {
                    Ok(batch_items) => results.extend(batch_items),
                    Err(e) => {
                        error!("Batch processing error: {:?}", e);
                        return Err(e); // Return the first processor error
                    }
                }
            }
            Err(panic_payload) => {
                let msg = panic_payload
                    .downcast_ref::<&'static str>()
                    .map(|s| *s)
                    .or_else(|| panic_payload.downcast_ref::<String>().map(|s| s.as_str()))
                    .unwrap_or("Unknown panic payload");
                error!("Thread panicked: {}", msg);
                return Err(ParallelError::ThreadError(format!("Thread panicked: {}", msg)).into());
            }
        }
    }

    Ok(results)
}

/// Helper struct wrapping a Rayon thread pool for executing parallel tasks.
pub struct ParallelExecutor {
    /// Rayon Thread pool instance
    pool: Arc<rayon::ThreadPool>,

    /// Configuration for tasks executed by this executor
    config: ParallelConfig,
}

impl ParallelExecutor {
    /// Create a new parallel executor with a configured Rayon thread pool.
    pub fn new(config: Option<ParallelConfig>) -> Result<Self, ParallelError> {
        let config = config.unwrap_or_default();

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(config.threads)
            .build()
            .map_err(|e| ParallelError::ThreadPoolBuildError(format!("{}", e)))?;

        Ok(ParallelExecutor {
            pool: Arc::new(pool),
            config,
        })
    }

    /// Execute an item-wise processing function in parallel using this executor's pool.
    ///
    /// Applies the `processor` function to each item in the `items` Vec using the
    /// configured Rayon thread pool and error handling strategy (`continue_on_error`).
    pub fn execute<T, U, F, E>(&self, items: Vec<T>, processor: F) -> Result<Vec<U>, E>
    where
        T: Send + Sync,
        U: Send,
        F: Fn(&T) -> Result<U, E> + Send + Sync,
        E: From<ParallelError> + Send + Debug,
    {
        self.pool
            .install(|| parallel_process(items, processor, Some(self.config.clone())))
    }
}
