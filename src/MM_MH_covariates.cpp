// src/MM_MH_covariates.cpp
#include <RcppArmadillo.h>
#include <cmath>
#include <random>
#include <algorithm>
#include <map>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// ----------------------------------------------------------------
// Helper functions
// ----------------------------------------------------------------

vec sample_dirichlet_cov(const vec& alpha, std::mt19937& rng) {
  int K = alpha.n_elem;
  vec samples(K);
  double sum_gamma = 0.0;
  for (int k = 0; k < K; ++k) {
    std::gamma_distribution<double> gamma_dist(alpha(k), 1.0);
    samples(k) = gamma_dist(rng);
    sum_gamma += samples(k);
  }
  if (sum_gamma > 0) samples = samples / sum_gamma;
  else samples.fill(1.0 / K);
  return samples;
}

double rinvgamma_cov(double shape, double scale, std::mt19937& rng) {
  std::gamma_distribution<double> gamma_dist(shape, 1.0 / scale);
  return 1.0 / gamma_dist(rng);
}

// ----------------------------------------------------------------
// Main MCMC sampler with adaptive proposal
// ----------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List run_MH_covariates(
    const arma::vec& y,
    const arma::ivec& category,
    const arma::mat& prior_means,
    int K,
    int n_iter = 10000,
    int burnin = 2000,
    double proposal_sd = 0.15,
    double alpha0 = 3.0,      // Stronger prior (was 0.001)
    double beta0 = 2.0,       // Stronger prior (was 0.001)
    double alpha_dirichlet = 5.0,
    int seed = 123,
    bool adaptive = true
) {
  // ---- Setup ----
  int N = y.n_elem;
  int C = prior_means.n_rows;
  int n_samples = n_iter - burnin;
  
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  std::normal_distribution<double> norm(0.0, 1.0);
  
  // ---- Initialize parameters ----
  mat mu = prior_means;
  mat sigma2 = ones(C, K) * 0.5;
  mat pi = ones(C, K) / K;
  
  // ---- Adaptive proposal ----
  mat proposal_sd_mat = ones(C, K) * proposal_sd;
  mat acceptance_rate(C, K, fill::zeros);
  int adaptation_window = 100;
  int total_adapt_iter = 0;
  double target_acceptance = 0.234;
  double adapt_c = 0.9;
  
  // ---- Initialize assignments ----
  imat z(N, K, fill::zeros);
  for (int i = 0; i < N; ++i) {
    int c = category(i);
    vec probs(K);
    for (int k = 0; k < K; ++k) {
      probs(k) = pi(c, k) * R::dnorm(y(i), mu(c, k), sqrt(sigma2(c, k)), 0);
    }
    probs /= sum(probs);
    double u = unif(rng);
    double cum = 0.0;
    for (int k = 0; k < K; ++k) {
      cum += probs(k);
      if (u <= cum) { z(i, k) = 1; break; }
    }
  }
  
  // ---- Storage ----
  cube mu_samples(C, K, n_samples);
  cube sigma2_samples(C, K, n_samples);
  cube pi_samples(C, K, n_samples);
  imat z_samples(N, n_samples, fill::zeros);
  
  // ---- MCMC loop ----
  int sample_idx = 0;
  
  for (int iter = 0; iter < n_iter; ++iter) {
    
    // 1. Update assignments (Gibbs)
    for (int i = 0; i < N; ++i) {
      int c = category(i);
      vec probs(K);
      for (int k = 0; k < K; ++k) {
        probs(k) = pi(c, k) * R::dnorm(y(i), mu(c, k), sqrt(sigma2(c, k)), 0);
      }
      probs /= sum(probs);
      double u = unif(rng);
      double cum = 0.0;
      z.row(i).zeros();
      for (int k = 0; k < K; ++k) {
        cum += probs(k);
        if (u <= cum) { z(i, k) = 1; break; }
      }
    }
    
    // 2. Update means (Metropolis-Hastings with adaptive proposal)
    for (int c = 0; c < C; ++c) {
      for (int k = 0; k < K; ++k) {
        // Count observations
        int n_ck = 0;
        for (int i = 0; i < N; ++i) {
          if (z(i, k) == 1 && category(i) == c) n_ck++;
        }
        if (n_ck == 0) {
          // No data: sample from prior
          mu(c, k) = prior_means(c, k) + norm(rng);
          continue;
        }
        
        double current_mu = mu(c, k);
        double current_sd = proposal_sd_mat(c, k);
        double proposed_mu = current_mu + current_sd * norm(rng);
        
        // Log-likelihood difference
        double log_lik_diff = 0.0;
        for (int i = 0; i < N; ++i) {
          if (z(i, k) == 1 && category(i) == c) {
            double yi = y(i);
            log_lik_diff += R::dnorm(yi, proposed_mu, sqrt(sigma2(c, k)), 1) -
              R::dnorm(yi, current_mu, sqrt(sigma2(c, k)), 1);
          }
        }
        
        // Prior difference
        double log_prior_diff = R::dnorm(proposed_mu, prior_means(c, k), 1.0, 1) -
          R::dnorm(current_mu, prior_means(c, k), 1.0, 1);
        
        double log_accept = log_lik_diff + log_prior_diff;
        
        if (log(unif(rng)) < log_accept) {
          mu(c, k) = proposed_mu;
          acceptance_rate(c, k) += 1.0;
        }
        total_adapt_iter++;
        
        // Adaptive tuning
        if (adaptive && (total_adapt_iter % adaptation_window == 0)) {
          double acc = acceptance_rate(c, k) / adaptation_window;
          double factor = exp(adapt_c * (acc - target_acceptance));
          proposal_sd_mat(c, k) = proposal_sd_mat(c, k) * std::max(0.01, std::min(10.0, factor));
          acceptance_rate(c, k) = 0.0;
        }
      }
    }
    
    // 3. Update variances (Gibbs) - FIXED with stronger prior to prevent Inf
    for (int c = 0; c < C; ++c) {
      for (int k = 0; k < K; ++k) {
        int n_ck = 0;
        double ss = 0.0;
        for (int i = 0; i < N; ++i) {
          if (z(i, k) == 1 && category(i) == c) {
            double yi = y(i);
            ss += (yi - mu(c, k)) * (yi - mu(c, k));
            n_ck++;
          }
        }
        
        // Use the stronger prior parameters passed from R
        if (n_ck > 0) {
          double shape = alpha0 + 0.5 * n_ck;
          double scale = beta0 + 0.5 * ss;
          sigma2(c, k) = rinvgamma_cov(shape, scale, rng);
        } else {
          // Empty component: sample from prior with strong regularization
          sigma2(c, k) = rinvgamma_cov(alpha0, beta0, rng);
          // Reinitialize mean to prior mean to help component recover
          mu(c, k) = prior_means(c, k) + 0.1 * norm(rng);
        }
      }
    }
    
    // 4. Update weights (Gibbs)
    for (int c = 0; c < C; ++c) {
      vec n_c_k(K, fill::zeros);
      for (int i = 0; i < N; ++i) {
        if (category(i) == c) {
          for (int k = 0; k < K; ++k) {
            if (z(i, k) == 1) { n_c_k(k) += 1.0; break; }
          }
        }
      }
      vec alpha_vec = ones(K) * alpha_dirichlet + n_c_k;
      pi.row(c) = sample_dirichlet_cov(alpha_vec, rng).t();
    }
    
    // 5. Fix label switching (ordering constraint)
    for (int c = 0; c < C; ++c) {
      // Get means for this category and sort them
      vec mu_c = mu.row(c).t();
      uvec order = sort_index(mu_c);
      
      // Reorder mu
      mu.row(c) = mu_c(order).t();
      
      // Reorder sigma2: convert to vector, reorder, put back
      vec sigma2_c = sigma2.row(c).t();
      sigma2.row(c) = sigma2_c(order).t();
      
      // Reorder pi: convert to vector, reorder, put back
      vec pi_c = pi.row(c).t();
      pi.row(c) = pi_c(order).t();
      
      // Update assignments for this category
      for (int i = 0; i < N; ++i) {
        if (category(i) == c) {
          for (int k = 0; k < K; ++k) {
            if (z(i, k) == 1) {
              z(i, k) = 0;
              z(i, order(k)) = 1;
              break;
            }
          }
        }
      }
    }
    
    // Store samples after burn-in
    if (iter >= burnin) {
      mu_samples.slice(sample_idx) = mu;
      sigma2_samples.slice(sample_idx) = sigma2;
      pi_samples.slice(sample_idx) = pi;
      for (int i = 0; i < N; ++i) {
        for (int k = 0; k < K; ++k) {
          if (z(i, k) == 1) {
            z_samples(i, sample_idx) = k + 1;
            break;
          }
        }
      }
      sample_idx++;
    }
  }
  
  // Return final proposal SD for diagnostics
  mat final_proposal_sd = mean(proposal_sd_mat);
  
  return Rcpp::List::create(
    Rcpp::Named("mu") = mu_samples,
    Rcpp::Named("sigma2") = sigma2_samples,
    Rcpp::Named("pi") = pi_samples,
    Rcpp::Named("z") = z_samples,
    Rcpp::Named("n_iter") = n_iter,
    Rcpp::Named("burnin") = burnin,
    Rcpp::Named("K") = K,
    Rcpp::Named("C") = C,
    Rcpp::Named("final_proposal_sd") = final_proposal_sd
  );
}