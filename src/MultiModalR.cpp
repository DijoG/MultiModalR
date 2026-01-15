// File: src/MultiModalR.cpp
#include <RcppArmadillo.h>
#include <cmath>
#include <random>
#include <algorithm>
#include <map>

// [[Rcpp::depends(RcppArmadillo)]]

// Helper function: sample from Dirichlet distribution
arma::vec rdirichlet_cpp(const arma::vec& alpha, std::mt19937& rng) {
  int K = alpha.n_elem;
  arma::vec draws(K);
  double sum_draws = 0.0;
  
  for(int k = 0; k < K; ++k) {
    std::gamma_distribution<double> gamma_dist(alpha(k), 1.0);
    draws(k) = gamma_dist(rng);
    sum_draws += draws(k);
  }
  
  return draws / sum_draws;
}

// Helper function: K-means++ initialization
arma::vec kmeans_plusplus(const arma::vec& y, int K, std::mt19937& rng) {
  int n = y.n_elem;
  arma::vec centers(K);
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  
  // First center: random point
  std::uniform_int_distribution<int> runif_int(0, n-1);
  centers(0) = y(runif_int(rng));
  
  for(int k = 1; k < K; ++k) {
    arma::vec distances(n);
    
    for(int i = 0; i < n; ++i) {
      double min_dist = INFINITY;
      for(int j = 0; j < k; ++j) {
        double dist = fabs(y(i) - centers(j));
        if(dist < min_dist) min_dist = dist;
      }
      distances(i) = min_dist * min_dist;
    }
    
    // Normalize to probabilities
    distances /= arma::sum(distances);
    
    // Cumulative sum
    arma::vec cumsum = arma::cumsum(distances);
    double u = unif(rng);
    
    int selected_idx = 0;
    for(int i = 0; i < n; ++i) {
      if(u <= cumsum(i)) {
        selected_idx = i;
        break;
      }
    }
    
    centers(k) = y(selected_idx);
  }
  
  return centers;
}

// MH-within-Gibbs implementation 
// Gaussian mixture MCMC that prevents empty groups
// [[Rcpp::export]]
Rcpp::List MM_MH_cpp(const arma::vec& y, 
                     const arma::vec& prior_means,
                     int n_iter = 1000, 
                     int burnin = 500,
                     double proposal_sd = 0.15,  // FOR MEAN UPDATES ONLY
                     int seed = 123) {
  
  int n = y.n_elem;
  int n_components = prior_means.n_elem;
  int n_samples = n_iter - burnin;
  
  // Random number generator
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  std::normal_distribution<double> norm(0.0, 1.0);
  std::uniform_int_distribution<int> runif_int(1, n_components);
  
  // Initialize storage for posterior probabilities
  arma::mat z_probs(n, n_components, arma::fill::zeros);
  
  // Initialize parameters
  arma::vec means = prior_means;
  arma::vec sds = arma::ones(n_components) * arma::stddev(y);
  arma::vec weights = arma::ones(n_components) / n_components;
  
  // Initialize z assignments - ensure all components have observations
  arma::ivec z_current(n);
  
  // Assign at least floor(n/n_components) to each component
  int min_per_component = std::max(1, n / n_components);
  int idx = 0;
  
  for(int k = 0; k < n_components; ++k) {
    for(int j = 0; j < min_per_component; ++j) {
      if(idx < n) {
        z_current(idx) = k + 1;
        idx++;
      }
    }
  }
  
  // Fill remaining with random assignments
  for(; idx < n; ++idx) {
    z_current(idx) = runif_int(rng);
  }
  
  // Shuffle to remove ordering bias
  for(int i = n - 1; i > 0; --i) {
    int j = std::rand() % (i + 1);
    std::swap(z_current(i), z_current(j));
  }
  
  // Acceptance tracking for mean updates
  arma::vec mean_acceptance(n_components, arma::fill::zeros);
  int total_mean_updates = 0;
  
  // MCMC loop
  for(int iter = 0; iter < n_iter; ++iter) {
    
    // ================= UPDATE COMPONENT MEANS =================
    for(int k = 0; k < n_components; ++k) {
      double mean_current = means(k);
      
      // PROPOSAL: use proposal_sd here
      double mean_proposed = mean_current + proposal_sd * norm(rng);
      
      double log_ratio = 0.0;
      int count = 0;
      
      for(int i = 0; i < n; ++i) {
        if(z_current(i) == (k + 1)) {
          double diff_prop = y(i) - mean_proposed;
          double diff_curr = y(i) - mean_current;
          log_ratio += (diff_curr * diff_curr - diff_prop * diff_prop) / (2.0 * sds(k) * sds(k));
          count++;
        }
      }
      
      // Prior: normal centered at prior_means[k]
      double prior_var = 1.0;
      log_ratio += (pow(mean_current - prior_means(k), 2) - 
        pow(mean_proposed - prior_means(k), 2)) / (2.0 * prior_var);
      
      // Accept/reject
      if(count > 0 && log(unif(rng)) < log_ratio) {
        means(k) = mean_proposed;
        mean_acceptance(k) += 1.0;
      }
      total_mean_updates++;
    }
    
    // ================= UPDATE GROUP ASSIGNMENTS =================
    // NOTE: Uses Gibbs sampling (categorical), NOT Metropolis with proposal_sd
    for(int i = 0; i < n; ++i) {
      // Calculate current counts (excluding current observation)
      arma::vec counts_excl_i(n_components, arma::fill::zeros);
      for(int j = 0; j < n; ++j) {
        if(j != i) {
          counts_excl_i(z_current(j) - 1) += 1.0;
        }
      }
      
      // Calculate log probabilities
      double max_log = -1e100;
      arma::vec log_probs(n_components);
      
      for(int k = 0; k < n_components; ++k) {
        double z_score = (y(i) - means(k)) / sds(k);
        
        // Log-likelihood term
        log_probs(k) = log(weights(k)) - 0.5 * z_score * z_score - log(sds(k));
        
        // Repellant term to prevent collapse
        double repellant_strength = 0.01;
        log_probs(k) -= repellant_strength * counts_excl_i(k);
        
        if(log_probs(k) > max_log) max_log = log_probs(k);
      }
      
      // Convert to probabilities
      double sum_exp = 0.0;
      for(int k = 0; k < n_components; ++k) {
        sum_exp += exp(log_probs(k) - max_log);
      }
      
      double inv_sum = 1.0 / sum_exp;
      arma::vec probs(n_components);
      
      for(int k = 0; k < n_components; ++k) {
        probs(k) = exp(log_probs(k) - max_log) * inv_sum;
      }
      
      // Sample from categorical distribution (Gibbs step)
      double u = unif(rng);
      double cum_prob = 0.0;
      int new_z = 1;
      
      for(int k = 0; k < n_components; ++k) {
        cum_prob += probs(k);
        if(u <= cum_prob) {
          new_z = k + 1;
          break;
        }
      }
      
      z_current(i) = new_z;
      
      // Store probabilities after burnin
      if(iter >= burnin) {
        for(int k = 0; k < n_components; ++k) {
          z_probs(i, k) += probs(k);
        }
      }
    }
    
    // Ensure no empty groups
    arma::vec counts(n_components, arma::fill::zeros);
    for(int i = 0; i < n; ++i) {
      counts(z_current(i) - 1) += 1.0;
    }
    
    for(int k = 0; k < n_components; ++k) {
      if(counts(k) == 0) {
        // Find largest group and move one observation
        int max_k = 0;
        int max_count = 0;
        for(int k2 = 0; k2 < n_components; ++k2) {
          if(counts(k2) > max_count) {
            max_count = counts(k2);
            max_k = k2;
          }
        }
        
        if(max_count > 1) {
          for(int i = 0; i < n; ++i) {
            if(z_current(i) == (max_k + 1)) {
              z_current(i) = k + 1;
              break;
            }
          }
        }
      }
    }
    
    // Update weights
    counts.zeros();
    for(int i = 0; i < n; ++i) {
      counts(z_current(i) - 1) += 1.0;
    }
    
    // Dirichlet posterior
    double dirichlet_alpha = 5.0;
    for(int k = 0; k < n_components; ++k) {
      weights(k) = (counts(k) + dirichlet_alpha) / (n + n_components * dirichlet_alpha);
    }
    
    // Update standard deviations occasionally
    if(iter % 20 == 0) {
      for(int k = 0; k < n_components; ++k) {
        double sum_sq = 0.0;
        int count_k = 0;
        
        for(int i = 0; i < n; ++i) {
          if(z_current(i) == (k + 1)) {
            double diff = y(i) - means(k);
            sum_sq += diff * diff;
            count_k++;
          }
        }
        
        if(count_k > 2) {
          // Inverse Gamma update
          double prior_shape = 3.0;
          double prior_rate = 1.0;
          double post_shape = prior_shape + count_k / 2.0;
          double post_rate = prior_rate + sum_sq / 2.0;
          
          sds(k) = sqrt(post_rate / (post_shape - 1));
          sds(k) = std::max(sds(k), 0.2);
        }
      }
    }
  }
  
  // Normalize probabilities
  z_probs = z_probs / n_samples;
  mean_acceptance = mean_acceptance / (total_mean_updates / n_components);
  
  // Calculate assignments
  arma::ivec assigned_group(n);
  for(int i = 0; i < n; ++i) {
    int best_comp = 1;  // ← START AT 1, NOT 0!
    double best_prob = z_probs(i, 0);  // Initialize with first component
    
    for(int k = 1; k < n_components; ++k) {  // Start from k=1
      if(z_probs(i, k) > best_prob) {
        best_prob = z_probs(i, k);
        best_comp = k + 1;
      }
    }
    
    // Additional safety: ensure best_comp is valid
    if(best_comp < 1 || best_comp > n_components) {
      best_comp = 1;  // Default to first component
    }
    
    assigned_group(i) = best_comp;
  }
  
  // Calculate group statistics
  arma::vec min_assigned(n), max_assigned(n), mean_assigned(n), mode_assigned(n);
  
  // Ensure all assigned_group values are valid
  for(int i = 0; i < n; ++i) {
    if(assigned_group(i) < 1 || assigned_group(i) > n_components) {
      assigned_group(i) = 1;
    }
  }
  
  for(int k = 0; k < n_components; ++k) {
    std::vector<double> group_values;
    
    for(int i = 0; i < n; ++i) {
      if(assigned_group(i) == (k + 1)) {
        group_values.push_back(y(i));
      }
    }
    
    if(!group_values.empty()) {
      double min_val = *std::min_element(group_values.begin(), group_values.end());
      double max_val = *std::max_element(group_values.begin(), group_values.end());
      
      double sum_val = 0.0;
      for(double val : group_values) sum_val += val;
      double mean_val = sum_val / group_values.size();
      
      // Simple mode estimation
      std::map<int, int> hist;
      int n_bins = std::min(20, static_cast<int>(group_values.size() / 5));
      if(n_bins < 3) n_bins = 3;
      
      double bin_width = (max_val - min_val) / n_bins;
      if(bin_width < 1e-10) bin_width = 0.1;
      
      for(double val : group_values) {
        int bin = static_cast<int>((val - min_val) / bin_width);
        if(bin >= 0 && bin < n_bins) hist[bin]++;
      }
      
      int mode_bin = 0;
      int max_count = 0;
      for(const auto& pair : hist) {
        if(pair.second > max_count) {
          max_count = pair.second;
          mode_bin = pair.first;
        }
      }
      
      double mode_val = min_val + (mode_bin + 0.5) * bin_width;
      
      for(int i = 0; i < n; ++i) {
        if(assigned_group(i) == (k + 1)) {
          min_assigned(i) = min_val;
          max_assigned(i) = max_val;
          mean_assigned(i) = mean_val;
          mode_assigned(i) = mode_val;
        }
      }
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("y") = y,
    Rcpp::Named("z_probs") = z_probs,
    Rcpp::Named("assigned_group") = assigned_group,
    Rcpp::Named("min_assigned") = min_assigned,
    Rcpp::Named("max_assigned") = max_assigned,
    Rcpp::Named("mean_assigned") = mean_assigned,
    Rcpp::Named("mode_assigned") = mode_assigned,
    Rcpp::Named("mean_acceptance") = mean_acceptance,
    Rcpp::Named("n_components") = n_components,
    Rcpp::Named("n_iter") = n_iter,
    Rcpp::Named("burnin") = burnin
  );
}

// Full Gibbs sampler for Gaussian mixture model with conjugate priors
// [[Rcpp::export]]
Rcpp::List MM_FullGibbs_cpp(const arma::vec& y,
                            const arma::vec& prior_means,
                            int n_iter = 1000,
                            int burnin = 500,
                            int seed = 123) {
  
  int n = y.n_elem;
  int K = prior_means.n_elem;
  int n_samples = n_iter - burnin;
  
  // Random number generator
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  std::normal_distribution<double> norm(0.0, 1.0);
  
  // ===================== PRIOR HYPERPARAMETERS =====================
  // Set reasonable defaults based on data
  double global_mean = arma::mean(y);
  double data_variance = arma::var(y);
  
  // For means: Normal(prior_mean, prior_tau2) - diffuse prior
  double prior_mean = global_mean;
  double prior_tau2 = 100.0 * data_variance;  // Diffuse
  
  // For variances: InverseGamma(alpha, beta) - weakly informative
  double alpha = 3.0;                         // Shape parameter
  double beta = data_variance;                // Rate parameter (scale)
  
  // For weights: Dirichlet(delta, ..., delta) - symmetric
  double delta = 5.0;                         // Concentration parameter
  
  // ===================== INITIALIZATION =====================
  
  // Initialize storage for posterior probabilities
  arma::mat z_probs(n, K, arma::fill::zeros);
  
  // Initialize parameters with K-means++ for better starting values
  arma::vec means = kmeans_plusplus(y, K, rng);
  arma::vec sds = arma::ones(K) * sqrt(data_variance);
  arma::vec weights = arma::ones(K) / K;
  
  // Initialize cluster assignments based on nearest mean
  arma::ivec z_current(n);
  
  for(int i = 0; i < n; ++i) {
    arma::vec dists(K);
    for(int k = 0; k < K; ++k) {
      dists(k) = fabs(y(i) - means(k));
    }
    z_current(i) = dists.index_min() + 1;  // 1-indexed
  }
  
  // Ensure no empty clusters initially
  arma::vec counts(K);
  for(int k = 0; k < K; ++k) {
    counts(k) = arma::sum(z_current == (k + 1));
    if(counts(k) == 0) {
      // Find largest cluster and split
      int max_k = counts.index_max();
      arma::uvec max_indices = arma::find(z_current == (max_k + 1));
      if(max_indices.n_elem > 1) {
        z_current(max_indices(0)) = k + 1;
      }
    }
  }
  
  // ===================== MCMC LOOP =====================
  
  for(int iter = 0; iter < n_iter; ++iter) {
    
    // ===== 1. UPDATE CLUSTER ASSIGNMENTS (Gibbs) =====
    for(int i = 0; i < n; ++i) {
      arma::vec log_probs(K);
      
      // Calculate log probabilities for each cluster
      for(int k = 0; k < K; ++k) {
        // Likelihood term (Gaussian)
        double z_score = (y(i) - means(k)) / sds(k);
        double log_lik = -0.5 * z_score * z_score - log(sds(k));
        
        // Prior term (log of weight)
        double log_prior = log(weights(k));
        
        // Full conditional (up to normalization)
        log_probs(k) = log_lik + log_prior;
      }
      
      // Numerical stability: subtract max
      double max_log = log_probs.max();
      log_probs -= max_log;
      
      // Convert to probabilities
      arma::vec probs = arma::exp(log_probs);
      probs /= arma::sum(probs);
      
      // Sample from categorical distribution
      double u = unif(rng);
      double cum_prob = 0.0;
      int new_z = 1;  // Default to cluster 1
      
      for(int k = 0; k < K; ++k) {
        cum_prob += probs(k);
        if(u <= cum_prob) {
          new_z = k + 1;
          break;
        }
      }
      
      z_current(i) = new_z;
      
      // Store probabilities after burnin
      if(iter >= burnin) {
        for(int k = 0; k < K; ++k) {
          z_probs(i, k) += probs(k);
        }
      }
    }
    
    // ===== 2. UPDATE CLUSTER COUNTS =====
    counts.zeros();
    for(int i = 0; i < n; ++i) {
      counts(z_current(i) - 1) += 1.0;
    }
    
    // ===== 3. UPDATE MEANS (Gibbs - Normal conjugate) =====
    for(int k = 0; k < K; ++k) {
      arma::uvec cluster_indices = arma::find(z_current == (k + 1));
      int n_k = cluster_indices.n_elem;
      
      if(n_k > 0) {
        // Data statistics
        double ybar_k = arma::mean(y.elem(cluster_indices));
        double sigma2_k = sds(k) * sds(k);
        
        // Prior precision and variance
        double prior_precision = 1.0 / prior_tau2;
        double data_precision = n_k / sigma2_k;
        
        // Posterior parameters
        double post_precision = prior_precision + data_precision;
        double post_variance = 1.0 / post_precision;
        double post_mean = (prior_precision * prior_mean + 
                            data_precision * ybar_k) / post_precision;
        
        // Sample from posterior
        means(k) = post_mean + sqrt(post_variance) * norm(rng);
      } else {
        // Sample from prior if cluster is empty
        means(k) = prior_mean + sqrt(prior_tau2) * norm(rng);
      }
    }
    
    // ===== 4. UPDATE VARIANCES (Gibbs - Inverse Gamma conjugate) =====
    for(int k = 0; k < K; ++k) {
      arma::uvec cluster_indices = arma::find(z_current == (k + 1));
      int n_k = cluster_indices.n_elem;
      
      if(n_k > 0) {
        // Calculate sum of squares
        double sum_sq = 0.0;
        for(int idx : cluster_indices) {
          double diff = y(idx) - means(k);
          sum_sq += diff * diff;
        }
        
        // Posterior parameters for Inverse Gamma
        double post_alpha = alpha + n_k / 2.0;
        double post_beta = beta + sum_sq / 2.0;
        
        // Sample from Inverse Gamma via Gamma
        std::gamma_distribution<double> gamma_dist(post_alpha, 1.0/post_beta);
        double inv_sigma2 = gamma_dist(rng);
        sds(k) = sqrt(1.0 / inv_sigma2);
        
        // Ensure minimum variance for numerical stability
        sds(k) = std::max(sds(k), 0.1);
      } else {
        // Sample from prior
        std::gamma_distribution<double> gamma_dist(alpha, 1.0/beta);
        double inv_sigma2 = gamma_dist(rng);
        sds(k) = sqrt(1.0 / inv_sigma2);
      }
    }
    
    // ===== 5. UPDATE WEIGHTS (Gibbs - Dirichlet conjugate) =====
    arma::vec dirichlet_params = counts + delta;
    weights = rdirichlet_cpp(dirichlet_params, rng);
    
    // Optional: Progress reporting
    if((iter + 1) % 1000 == 0) {
      Rcpp::Rcout << "Iteration " << (iter + 1) << "/" << n_iter;
      Rcpp::Rcout << " | Non-empty clusters: " << arma::sum(counts > 0);
      Rcpp::Rcout << std::endl;
    }
  }
  
  // ===================== POST-PROCESSING =====================
  
  // Normalize probabilities
  z_probs = z_probs / n_samples;
  
  // Calculate assignments and confidence
  arma::ivec assigned_group(n);
  arma::vec assignment_confidence(n);
  
  for(int i = 0; i < n; ++i) {
    double max_prob = z_probs(i, 0);
    int best_k = 0;
    
    for(int k = 1; k < K; ++k) {
      if(z_probs(i, k) > max_prob) {
        max_prob = z_probs(i, k);
        best_k = k;
      }
    }
    
    assigned_group(i) = best_k + 1;  // 1-indexed
    assignment_confidence(i) = max_prob;
  }
  
  // Calculate group statistics
  arma::vec min_assigned(n), max_assigned(n), mean_assigned(n), mode_assigned(n);
  
  for(int k = 0; k < K; ++k) {
    arma::uvec group_indices = arma::find(assigned_group == (k + 1));
    
    if(group_indices.n_elem > 0) {
      arma::vec group_vals = y.elem(group_indices);
      
      double min_val = arma::min(group_vals);
      double max_val = arma::max(group_vals);
      double mean_val = arma::mean(group_vals);
      
      // Simple mode estimation
      int n_bins = std::min(20, static_cast<int>(group_vals.n_elem / 5));
      n_bins = std::max(n_bins, 3);
      
      arma::vec hist_counts(n_bins, arma::fill::zeros);
      double bin_width = (max_val - min_val) / n_bins;
      
      if(bin_width < 1e-10) {
        bin_width = 0.1;
      }
      
      for(int idx = 0; idx < group_vals.n_elem; ++idx) {
        int bin_idx = static_cast<int>((group_vals(idx) - min_val) / bin_width);
        bin_idx = std::min(std::max(bin_idx, 0), n_bins - 1);
        hist_counts(bin_idx)++;
      }
      
      int mode_bin = hist_counts.index_max();
      double mode_val = min_val + (mode_bin + 0.5) * bin_width;
      
      // Assign to all members of this group
      for(int idx = 0; idx < group_indices.n_elem; ++idx) {
        int i = group_indices(idx);
        min_assigned(i) = min_val;
        max_assigned(i) = max_val;
        mean_assigned(i) = mean_val;
        mode_assigned(i) = mode_val;
      }
    }
  }
  
  // Calculate posterior means for parameters
  arma::vec posterior_means = means;
  arma::vec posterior_sds = sds;
  arma::vec posterior_weights = weights;
  
  // ===================== RETURN RESULTS =====================
  
  return Rcpp::List::create(
    Rcpp::Named("y") = y,
    Rcpp::Named("z_probs") = z_probs,
    Rcpp::Named("assigned_group") = assigned_group,
    Rcpp::Named("assignment_confidence") = assignment_confidence,
    Rcpp::Named("min_assigned") = min_assigned,
    Rcpp::Named("max_assigned") = max_assigned,
    Rcpp::Named("mean_assigned") = mean_assigned,
    Rcpp::Named("mode_assigned") = mode_assigned,
    Rcpp::Named("posterior_means") = posterior_means,
    Rcpp::Named("posterior_sds") = posterior_sds,
    Rcpp::Named("posterior_weights") = posterior_weights,
    Rcpp::Named("n_components") = K,
    Rcpp::Named("n_iter") = n_iter,
    Rcpp::Named("burnin") = burnin,
    Rcpp::Named("prior_hyperparameters") = Rcpp::List::create(
      Rcpp::Named("prior_mean") = prior_mean,
      Rcpp::Named("prior_tau2") = prior_tau2,
      Rcpp::Named("alpha") = alpha,
      Rcpp::Named("beta") = beta,
      Rcpp::Named("delta") = delta
    )
  );
}