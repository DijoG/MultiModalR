// src/MM_MH_covariates.cpp
#include <RcppArmadillo.h>
#include <cmath>
#include <random>
#include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// ----------------------------------------------------------------
// Helper functions
// ----------------------------------------------------------------

vec rdirichlet(const vec& alpha, std::mt19937& rng) {
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

double rinvgamma(double shape, double scale, std::mt19937& rng) {
  std::gamma_distribution<double> gamma_dist(shape, 1.0 / scale);
  return 1.0 / gamma_dist(rng);
}

// ----------------------------------------------------------------
// Main MCMC sampler with PER-CATEGORY K
// ----------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List run_MH_covariates(
    const arma::vec& y,
    const arma::ivec& category,
    const arma::ivec& K_per_category,   // <-- NEW
    const arma::mat& prior_means,
    int maxK,                           // <-- NEW (was K)
    int n_iter = 10000,
    int burnin = 2000,
    double proposal_sd = 0.15,
    double alpha0 = 3.0,
    double beta0 = 2.0,
    double alpha_dirichlet = 5.0,
    int seed = 123,
    bool adaptive = true
) {
  // ---- Setup ----
  int N = y.n_elem;
  int C = K_per_category.n_elem;
  int n_samples = n_iter - burnin;
  
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> unif(0.0, 1.0);
  std::normal_distribution<double> norm(0.0, 1.0);
  
  // ---- Initialize parameters ----
  // mu: C x maxK (but only K_per_category[c] are used for each category)
  mat mu = prior_means;
  mat sigma2 = ones(C, maxK) * 0.5;
  mat pi = ones(C, maxK) / maxK;
  
  // ---- Initialize assignments ----
  imat z(N, maxK, fill::zeros);
  for (int i = 0; i < N; ++i) {
    int c = category(i);
    int K_c = K_per_category(c);
    vec probs(maxK, fill::zeros);
    double sum_probs = 0.0;
    for (int k = 0; k < K_c; ++k) {
      probs(k) = pi(c, k) * R::dnorm(y(i), mu(c, k), sqrt(sigma2(c, k)), 0);
      sum_probs += probs(k);
    }
    if (sum_probs > 0) {
      probs /= sum_probs;
    } else {
      probs.fill(1.0 / K_c);
    }
    double u = unif(rng);
    double cum = 0.0;
    for (int k = 0; k < K_c; ++k) {
      cum += probs(k);
      if (u <= cum) {
        z(i, k) = 1;
        break;
      }
    }
  }
  
  // ---- Storage ----
  cube mu_samples(C, maxK, n_samples);
  cube sigma2_samples(C, maxK, n_samples);
  cube pi_samples(C, maxK, n_samples);
  imat z_samples(N, n_samples, fill::zeros);
  
  // ---- Adaptive proposal ----
  mat proposal_sd_mat = ones(C, maxK) * proposal_sd;
  mat acceptance_rate(C, maxK, fill::zeros);
  int adaptation_window = 100;
  int total_adapt_iter = 0;
  double target_acceptance = 0.234;
  double adapt_c = 0.9;
  
  // ---- MCMC loop ----
  int sample_idx = 0;
  
  for (int iter = 0; iter < n_iter; ++iter) {
    
    // 1. Update assignments
    for (int i = 0; i < N; ++i) {
      int c = category(i);
      int K_c = K_per_category(c);
      vec probs(maxK, fill::zeros);
      double sum_probs = 0.0;
      for (int k = 0; k < K_c; ++k) {
        probs(k) = pi(c, k) * R::dnorm(y(i), mu(c, k), sqrt(sigma2(c, k)), 0);
        sum_probs += probs(k);
      }
      if (sum_probs > 0) {
        probs /= sum_probs;
      } else {
        probs.fill(1.0 / K_c);
      }
      double u = unif(rng);
      double cum = 0.0;
      z.row(i).zeros();
      for (int k = 0; k < K_c; ++k) {
        cum += probs(k);
        if (u <= cum) {
          z(i, k) = 1;
          break;
        }
      }
    }
    
    // 2. Update means (Metropolis-Hastings) - PER CATEGORY
    for (int c = 0; c < C; ++c) {
      int K_c = K_per_category(c);
      for (int k = 0; k < K_c; ++k) {
        // Count observations in category c assigned to component k
        int n_ck = 0;
        for (int i = 0; i < N; ++i) {
          if (z(i, k) == 1 && category(i) == c) n_ck++;
        }
        if (n_ck == 0) {
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
        
        if (adaptive && (total_adapt_iter % adaptation_window == 0)) {
          double acc = acceptance_rate(c, k) / adaptation_window;
          double factor = exp(adapt_c * (acc - target_acceptance));
          proposal_sd_mat(c, k) = proposal_sd_mat(c, k) * std::max(0.01, std::min(10.0, factor));
          acceptance_rate(c, k) = 0.0;
        }
      }
    }
    
    // 3. Update variances (Gibbs) - PER CATEGORY
    for (int c = 0; c < C; ++c) {
      int K_c = K_per_category(c);
      for (int k = 0; k < K_c; ++k) {
        int n_ck = 0;
        double ss = 0.0;
        for (int i = 0; i < N; ++i) {
          if (z(i, k) == 1 && category(i) == c) {
            double yi = y(i);
            ss += (yi - mu(c, k)) * (yi - mu(c, k));
            n_ck++;
          }
        }
        if (n_ck > 0) {
          double shape = alpha0 + 0.5 * n_ck;
          double scale = beta0 + 0.5 * ss;
          sigma2(c, k) = rinvgamma(shape, scale, rng);
        } else {
          sigma2(c, k) = rinvgamma(alpha0, beta0, rng);
        }
      }
    }
    
    // 4. Update weights (Gibbs) - PER CATEGORY
    for (int c = 0; c < C; ++c) {
      int K_c = K_per_category(c);
      vec n_c_k(K_c, fill::zeros);
      for (int i = 0; i < N; ++i) {
        if (category(i) == c) {
          for (int k = 0; k < K_c; ++k) {
            if (z(i, k) == 1) { n_c_k(k) += 1.0; break; }
          }
        }
      }
      vec alpha_vec = ones(K_c) * alpha_dirichlet + n_c_k;
      vec pi_c = rdirichlet(alpha_vec, rng);
      for (int k = 0; k < K_c; ++k) {
        pi(c, k) = pi_c(k);
      }
      // Set unused components to 0
      for (int k = K_c; k < maxK; ++k) {
        pi(c, k) = 0.0;
      }
    }
    
    // 5. Fix label switching (ordering constraint PER CATEGORY)
    for (int c = 0; c < C; ++c) {
      int K_c = K_per_category(c);
      if (K_c > 1) {
        // Get means for this category
        vec mu_c = mu.row(c).t();
        uvec order = sort_index(mu_c.subvec(0, K_c - 1));
        
        // Reorder mu, sigma2, pi for this category
        mat mu_c_new = mu.row(c);
        mat sigma2_c_new = sigma2.row(c);
        mat pi_c_new = pi.row(c);
        for (int k = 0; k < K_c; ++k) {
          mu(c, k) = mu_c_new(order(k));
          sigma2(c, k) = sigma2_c_new(order(k));
          pi(c, k) = pi_c_new(order(k));
        }
        
        // Update assignments for this category
        for (int i = 0; i < N; ++i) {
          if (category(i) == c) {
            for (int k = 0; k < K_c; ++k) {
              if (z(i, k) == 1) {
                z(i, k) = 0;
                z(i, order(k)) = 1;
                break;
              }
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
        for (int k = 0; k < maxK; ++k) {
          if (z(i, k) == 1) {
            z_samples(i, sample_idx) = k + 1;
            break;
          }
        }
      }
      sample_idx++;
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("mu") = mu_samples,
    Rcpp::Named("sigma2") = sigma2_samples,
    Rcpp::Named("pi") = pi_samples,
    Rcpp::Named("z") = z_samples,
    Rcpp::Named("n_iter") = n_iter,
    Rcpp::Named("burnin") = burnin,
    Rcpp::Named("maxK") = maxK,
    Rcpp::Named("C") = C,
    Rcpp::Named("K_per_category") = K_per_category
  );
}