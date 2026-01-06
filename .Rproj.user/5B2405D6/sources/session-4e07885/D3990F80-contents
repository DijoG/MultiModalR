// File: src/MultiModalR.cpp
#include <RcppArmadillo.h>
#include <cmath>
#include <random>
#include <algorithm>
#include <map>

// [[Rcpp::depends(RcppArmadillo)]]

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
    int best_comp = 0;
    double best_prob = 0.0;
    for(int k = 0; k < n_components; ++k) {
      if(z_probs(i, k) > best_prob) {
        best_prob = z_probs(i, k);
        best_comp = k + 1;
      }
    }
    assigned_group(i) = best_comp;
  }
  
  // Calculate group statistics
  arma::vec min_assigned(n), max_assigned(n), mean_assigned(n), mode_assigned(n);
  
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