use nalgebra::ComplexField;
use ndarray::{Array1, Array2};
use rand::prelude::*;
// Import random libraries with feature flag
#[cfg(feature = "random")]
use rand_distr::{Dirichlet, Distribution};
use std::collections::HashMap;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum BayesianError {
    #[error("Invalid dimensions: observed features {0} != signature features {1}")]
    DimensionMismatch(usize, usize),
    #[error("MCMC convergence failure after {0} iterations")]
    ConvergenceFailure(usize),
}

/// Results from strain abundance estimation
pub struct StrainAbundanceResult {
    /// Strain ID to (abundance, confidence interval) mapping
    pub abundances: HashMap<String, (f64, f64)>,
    /// Effective sample size from MCMC
    pub effective_sample_size: f64,
    /// Model fit quality measure
    pub goodness_of_fit: f64,
}

/// Bayesian mixture model for strain deconvolution
pub struct StrainMixtureModel {
    /// Number of strains in the model
    n_strains: usize,
    /// Number of features (k-mer counts, etc.)
    n_features: usize,
    /// Strain signature matrix (features x strains)
    signatures: Array2<f64>,
    /// Strain IDs
    strain_ids: Vec<String>,
    /// Prior distribution parameters for strain abundances
    abundance_prior: Vec<f64>,
    /// MCMC parameters
    mcmc_iterations: usize,
    mcmc_burnin: usize,
    mcmc_thin: usize,
    /// Random number generator
    rng: StdRng,
}

impl StrainMixtureModel {
    /// Create a new strain mixture model
    pub fn new(
        signatures: Array2<f64>,
        strain_ids: Vec<String>,
        abundance_prior: Option<Vec<f64>>,
        mcmc_iterations: Option<usize>,
        seed: Option<u64>,
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let n_strains = signatures.shape()[1];
        let n_features = signatures.shape()[0];

        // Validate inputs
        if strain_ids.len() != n_strains {
            return Err(Box::new(BayesianError::DimensionMismatch(
                strain_ids.len(),
                n_strains,
            )));
        }

        // Create default prior if not provided
        let prior = match abundance_prior {
            Some(p) => {
                if p.len() != n_strains {
                    return Err(Box::new(BayesianError::DimensionMismatch(
                        p.len(),
                        n_strains,
                    )));
                }
                p
            }
            None => vec![1.0; n_strains], // Uniform Dirichlet prior
        };

        // Set MCMC parameters with reasonable defaults
        let iterations = mcmc_iterations.unwrap_or(10000);

        // Initialize RNG with provided seed or random seed
        #[cfg(feature = "random")]
        let random_seed = seed.unwrap_or_else(|| rand::thread_rng().gen());

        #[cfg(not(feature = "random"))]
        let random_seed = seed.unwrap_or(0);

        let rng = StdRng::seed_from_u64(random_seed);

        Ok(Self {
            n_strains,
            n_features,
            signatures,
            strain_ids,
            abundance_prior: prior,
            mcmc_iterations: iterations,
            mcmc_burnin: iterations / 5, // 20% burnin
            mcmc_thin: 10,               // Keep every 10th sample
            rng,
        })
    }

    /// Estimate strain abundances from observed data using MCMC
    pub fn estimate_abundances(
        &mut self,
        observed: &Array1<f64>,
    ) -> Result<StrainAbundanceResult, Box<dyn std::error::Error>> {
        // Validate dimensions
        if observed.len() != self.n_features {
            return Err(Box::new(BayesianError::DimensionMismatch(
                observed.len(),
                self.n_features,
            )));
        }

        // Normalize observed data
        let total = observed.sum();
        let observed_norm = if total > 0.0 {
            observed / total
        } else {
            observed.clone()
        };

        // Initialize abundance vector (starting point for MCMC)
        let mut current_abundances =
            Array1::from_vec(vec![1.0 / self.n_strains as f64; self.n_strains]);

        // Storage for MCMC samples
        let samples_to_store = (self.mcmc_iterations - self.mcmc_burnin) / self.mcmc_thin;
        let mut abundance_samples = Vec::with_capacity(samples_to_store);
        let mut log_likelihood_samples = Vec::with_capacity(samples_to_store);

        // Run MCMC
        // Simplified MCMC implementation that doesn't require Dirichlet
        for iteration in 0..self.mcmc_iterations {
            // Propose new abundances using simple perturbation
            let mut proposal = Vec::with_capacity(self.n_strains);
            let perturbation_scale = 0.1;

            // Generate perturbed values
            for a in &current_abundances {
                let mut new_val = a + self.rng.gen_range(-perturbation_scale..perturbation_scale);
                if new_val < 0.0 {
                    new_val = 0.0;
                }
                proposal.push(new_val);
            }

            // Normalize to sum to 1
            let sum: f64 = proposal.iter().sum();
            if sum > 0.0 {
                for val in &mut proposal {
                    *val /= sum;
                }
            }

            let proposal_abundances = Array1::from_vec(proposal);

            // Calculate likelihood for current and proposed abundances
            let current_likelihood = self.calculate_likelihood(&observed_norm, &current_abundances);
            let proposal_likelihood =
                self.calculate_likelihood(&observed_norm, &proposal_abundances);

            // Accept/reject based on likelihood ratio (simplified Metropolis-Hastings)
            let log_acceptance_ratio = proposal_likelihood - current_likelihood;
            let random_log = (self.rng.random_range(0.0..1.0) + 1e-10).ln(); // Add small value to avoid ln(0)
            let accept = log_acceptance_ratio > 0.0 || random_log < log_acceptance_ratio;

            if accept {
                current_abundances = proposal_abundances;
            }

            // Store samples after burn-in, applying thinning
            if iteration >= self.mcmc_burnin && (iteration - self.mcmc_burnin) % self.mcmc_thin == 0
            {
                abundance_samples.push(current_abundances.to_vec());
                log_likelihood_samples.push(current_likelihood);
            }
        }

        // Process MCMC samples to get abundance estimates and confidence intervals
        let result = self.process_samples(abundance_samples, log_likelihood_samples);

        Ok(result)
    }

    /// Calculate log-likelihood of observed data given abundance parameters
    fn calculate_likelihood(&self, observed: &Array1<f64>, abundances: &Array1<f64>) -> f64 {
        // Generate expected profile by mixing strain signatures according to abundances
        let expected = self.signatures.dot(abundances);

        // Calculate negative log-likelihood (assuming Poisson model for counts)
        let mut log_likelihood = 0.0;

        for (obs, exp) in observed.iter().zip(expected.iter()) {
            if *exp > 0.0 {
                // Poisson log-likelihood (ignoring constants)
                log_likelihood += obs * exp.ln() - *exp;
            }
        }

        log_likelihood
    }

    /// Process MCMC samples to get abundance estimates and confidence intervals
    fn process_samples(
        &self,
        abundance_samples: Vec<Vec<f64>>,
        likelihood_samples: Vec<f64>,
    ) -> StrainAbundanceResult {
        let n_samples = abundance_samples.len();
        let mut abundances = HashMap::new();

        // Calculate mean and 95% confidence interval for each strain
        for i in 0..self.n_strains {
            let strain_id = &self.strain_ids[i];

            // Extract all samples for this strain
            let mut strain_samples: Vec<f64> =
                abundance_samples.iter().map(|sample| sample[i]).collect();

            // Sort for percentile calculation
            strain_samples.sort_by(|a, b| a.partial_cmp(b).unwrap());

            // Calculate mean abundance
            let mean_abundance = strain_samples.iter().sum::<f64>() / n_samples as f64;

            // Calculate confidence interval width (using 95% credible interval)
            let lower_idx = (0.025 * n_samples as f64) as usize;
            let upper_idx = (0.975 * n_samples as f64) as usize;
            let confidence_interval = strain_samples[upper_idx] - strain_samples[lower_idx];

            abundances.insert(strain_id.clone(), (mean_abundance, confidence_interval));
        }

        // Calculate effective sample size (simplified)
        let effective_sample_size = n_samples as f64;

        // Calculate goodness of fit (using mean log-likelihood)
        let goodness_of_fit =
            likelihood_samples.iter().sum::<f64>() / likelihood_samples.len() as f64;

        StrainAbundanceResult {
            abundances,
            effective_sample_size,
            goodness_of_fit,
        }
    }
}

/// Strain deconvolution algorithm for metagenomic samples
///
/// This struct provides methods for estimating the relative abundance
/// of different microbial strains in a mixed sample.
#[derive(Debug)]
pub struct StrainDeconvolution {
    /// Reference strain signatures
    pub reference_signatures: Vec<Array1<f64>>,

    /// Reference strain IDs
    pub reference_ids: Vec<String>,

    /// Minimum abundance threshold to report
    pub min_abundance: f64,

    /// Maximum iterations for optimization
    pub max_iterations: usize,
}

impl StrainDeconvolution {
    /// Create a new StrainDeconvolution instance
    pub fn new(
        reference_signatures: Vec<Array1<f64>>,
        reference_ids: Vec<String>,
        min_abundance: Option<f64>,
        max_iterations: Option<usize>,
    ) -> Result<Self, String> {
        // Validate inputs
        if reference_signatures.len() != reference_ids.len() {
            return Err(format!(
                "Number of signatures ({}) does not match number of IDs ({})",
                reference_signatures.len(),
                reference_ids.len()
            ));
        }

        if reference_signatures.is_empty() {
            return Err("No reference signatures provided".to_string());
        }

        Ok(Self {
            reference_signatures,
            reference_ids,
            min_abundance: min_abundance.unwrap_or(0.01), // Default 1%
            max_iterations: max_iterations.unwrap_or(1000),
        })
    }

    /// Estimate strain abundances in a sample using NNLS
    ///
    /// # Arguments
    ///
    /// * `sample_profile` - Feature vector from metagenomic sample
    ///
    /// # Returns
    ///
    /// HashMap mapping strain IDs to their estimated abundances
    pub fn estimate_abundances(&self, sample_profile: &Array1<f64>) -> HashMap<String, f64> {
        // This is a simplified implementation using basic non-negative least squares
        // In practice, would likely use NNLS from an optimized library

        // Build signature matrix from reference signatures
        let n_features = sample_profile.len();
        let n_strains = self.reference_signatures.len();

        let mut signature_matrix = Array2::<f64>::zeros((n_features, n_strains));
        for (i, sig) in self.reference_signatures.iter().enumerate() {
            let mut col = signature_matrix.column_mut(i);
            col.assign(sig);
        }

        // Initialize abundance vector to equal abundances
        let mut abundances = Array1::<f64>::ones(n_strains) / n_strains as f64;

        // Placeholder for actual optimization
        // In a real implementation, would use an NNLS solver
        for _ in 0..self.max_iterations {
            // Compute current prediction
            let prediction = signature_matrix.dot(&abundances);

            // Compute residual
            let residual = sample_profile - &prediction;

            // Compute gradient (simplified)
            let gradient = signature_matrix.t().dot(&residual);

            // Update abundances with a small step in gradient direction
            let step_size = 0.01;
            abundances = &abundances + step_size * gradient;

            // Project to non-negative values
            for a in abundances.iter_mut() {
                if *a < 0.0 {
                    *a = 0.0;
                }
            }

            // Normalize to sum to 1
            let sum = abundances.sum();
            if sum > 0.0 {
                abundances /= sum;
            }
        }

        // Convert to HashMap, filtering by minimum abundance
        let mut result = HashMap::new();
        for (i, strain_id) in self.reference_ids.iter().enumerate() {
            let abundance = abundances[i];
            if abundance >= self.min_abundance {
                result.insert(strain_id.clone(), abundance);
            }
        }

        result
    }
}
