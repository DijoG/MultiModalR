# R/fuss_COVARIATE_mcmc.R

#' Unified Bayesian Mixture Model with Category-Specific Parameters
#'
#' Fits a Bayesian Gaussian mixture model where each category can have
#' a different number of components (K_c), while sharing information
#' across categories for means and variances.
#'
#' @param data Data frame containing the variables.
#' @param varY Name of the continuous response variable.
#' @param varCLASS Name of the categorical grouping variable.
#' @param varID Optional name of the ID variable.
#' @param K Optional vector of K per category. If NULL, auto-detected.
#' @param n_iter Number of MCMC iterations (default: 10000).
#' @param burnin Number of burn-in iterations (default: 2000).
#' @param proposal_sd Proposal standard deviation (default: 0.15).
#' @param adaptive Logical; adapt proposal variance (default: TRUE).
#' @param alpha0 Shape parameter for Inverse Gamma prior (default: 3.0).
#' @param beta0 Scale parameter for Inverse Gamma prior (default: 2.0).
#' @param alpha_dirichlet Concentration parameter (default: 5.0).
#' @param method Bandwidth selection method (default: "sj-dpi").
#' @param sj_adjust Adjustment factor (default: 0.5).
#' @param within Merging radius (default: 1.0).
#' @param seed Random seed (default: 123).
#'
#' @return A list with components.
#' @export
fuss_COVARIATE_mcmc <- function(
    data,
    varY,
    varCLASS,
    varID = NULL,
    K = NULL,
    n_iter = 10000,
    burnin = 2000,
    proposal_sd = 0.15,
    adaptive = TRUE,
    alpha0 = 3.0,
    beta0 = 2.0,
    alpha_dirichlet = 5.0,
    method = "sj-dpi",
    sj_adjust = 0.5,
    within = 1.0,
    seed = 123
) {
  # ---- Input validation ----
  if (!is.data.frame(data)) stop("'data' must be a data frame")
  if (!varY %in% names(data)) stop("'varY' not found in data")
  if (!varCLASS %in% names(data)) stop("'varCLASS' not found in data")
  
  y = data[[varY]]
  category = data[[varCLASS]]
  if (!is.numeric(y)) stop("'varY' must be numeric")
  
  # Convert category to numeric (0-based for C++)
  cat_factor = as.factor(category)
  C = nlevels(cat_factor)
  cat_numeric = as.numeric(cat_factor) - 1
  
  # ---- Detect K per category ----
  if (is.null(K)) {
    K_per_category = numeric(C)
    for (c in 1:C) {
      y_c = y[cat_numeric == (c - 1)]
      if (length(y_c) > 5) {
        cat_modes_list = get_MODES_enhanced(y_c, adjust = sj_adjust, threshold = 1.0)
        modes_df = cat_modes_list[[method]]
        if (!is.null(modes_df) && nrow(modes_df) > 0) {
          grouped_modes = group_MODES_enhanced(modes_df, within = within)
          K_per_category[c] = max(2, length(grouped_modes$Est_Mode))
        } else {
          K_per_category[c] = 3
        }
      } else {
        K_per_category[c] = 2
      }
    }
    K = as.integer(K_per_category)
    message("Detected K per category: ", paste(K, collapse = ", "))
  } else {
    if (length(K) != C) {
      stop("'K' must be a vector of length C (number of categories = ", C, ")")
    }
    K = as.integer(K)
  }
  
  maxK = max(K)
  
  # ---- Build prior means matrix (C x maxK) ----
  prior_means = matrix(NA, C, maxK)
  for (c in 1:C) {
    y_c = y[cat_numeric == (c - 1)]
    K_c = K[c]
    if (length(y_c) > 5) {
      cat_modes_list = get_MODES_enhanced(y_c, adjust = sj_adjust, threshold = 1.0)
      modes_df = cat_modes_list[[method]]
      if (!is.null(modes_df) && nrow(modes_df) > 0) {
        grouped_modes = group_MODES_enhanced(modes_df, within = within)
        modes = grouped_modes$Est_Mode
      } else {
        modes = quantile(y_c, probs = seq(0.2, 0.8, length.out = maxK))
      }
    } else {
      modes = quantile(y_c, probs = seq(0.2, 0.8, length.out = maxK))
    }
    
    if (length(modes) >= maxK) {
      prior_means[c, ] = sort(modes[1:maxK])
    } else {
      prior_means[c, ] = quantile(y_c, probs = seq(0.1, 0.9, length.out = maxK))
    }
  }
  
  # ---- Call C++ sampler ----
  cpp_result = MultiModalR:::run_MH_covariates(
    y = y,
    category = cat_numeric,
    K_per_category = K,
    prior_means = prior_means,
    maxK = maxK,
    n_iter = n_iter,
    burnin = burnin,
    proposal_sd = proposal_sd,
    alpha0 = alpha0,
    beta0 = beta0,
    alpha_dirichlet = alpha_dirichlet,
    seed = seed,
    adaptive = adaptive
  )
  
  # ---- Post-process ----
  mu_samples = cpp_result$mu
  sigma2_samples = cpp_result$sigma2
  pi_samples = cpp_result$pi
  z_samples = cpp_result$z
  
  N = length(y)
  N_samples = dim(mu_samples)[3]
  maxK = dim(mu_samples)[2]
  
  # Posterior means
  mu_mean = apply(mu_samples, c(1, 2), mean)
  sigma2_mean = apply(sigma2_samples, c(1, 2), mean)
  pi_mean = apply(pi_samples, c(1, 2), mean)
  
  # Most likely assignment (mode across samples)
  assignments = apply(z_samples, 1, function(x) {
    tab = table(x)
    as.numeric(names(tab)[which.max(tab)])
  })
  
  # ---- Compute per-observation probability matrix ----
  prob_matrix = matrix(0, N, maxK)
  
  for (s in 1:N_samples) {
    mu_s = mu_samples[, , s]
    sigma2_s = sigma2_samples[, , s]
    pi_s = pi_samples[, , s]
    
    if (is.vector(mu_s)) {
      mu_s = matrix(mu_s, nrow = C, ncol = maxK)
    }
    if (is.vector(sigma2_s)) {
      sigma2_s = matrix(sigma2_s, nrow = C, ncol = maxK)
    }
    if (is.vector(pi_s)) {
      pi_s = matrix(pi_s, nrow = C, ncol = maxK)
    }
    
    for (i in 1:N) {
      c_idx = cat_numeric[i] + 1
      # Only consider components up to K_per_category[c_idx]
      K_c = K[c_idx]
      for (k in 1:K_c) {
        prob_matrix[i, k] = prob_matrix[i, k] + 
          pi_s[c_idx, k] * dnorm(y[i], mu_s[c_idx, k], sqrt(sigma2_s[c_idx, k]))
      }
    }
  }
  
  prob_matrix = prob_matrix / N_samples
  row_sums = rowSums(prob_matrix)
  if (any(row_sums == 0)) {
    prob_matrix = prob_matrix + 1e-10
    row_sums = rowSums(prob_matrix)
  }
  prob_matrix = prob_matrix / row_sums
  
  colnames(prob_matrix) = paste0("Group_", 1:maxK)
  prob_df = as.data.frame(prob_matrix)
  
  # ---- Compute group statistics ----
  min_assigned = numeric(N)
  max_assigned = numeric(N)
  mean_assigned = numeric(N)
  mode_assigned = numeric(N)
  
  for (k in 1:maxK) {
    idx = assignments == k
    n_k = sum(idx)
    if (n_k > 0) {
      y_k = y[idx]
      min_assigned[idx] = min(y_k)
      max_assigned[idx] = max(y_k)
      mean_assigned[idx] = mean(y_k)
      
      if (n_k >= 3) {
        tryCatch({
          dens = density(y_k, n = 128)
          mode_assigned[idx] = dens$x[which.max(dens$y)]
        }, error = function(e) {
          mode_assigned[idx] = mean(y_k)
        })
      } else {
        mode_assigned[idx] = mean(y_k)
      }
    }
  }
  
  # ---- Build output ----
  out = list(
    y = y,
    ID = if (!is.null(varID)) data[[varID]] else 1:N,
    Main_Class = data[[varCLASS]],
    prob_matrix = prob_df,
    Assigned_Group = assignments,
    Min_Assigned = min_assigned,
    Max_Assigned = max_assigned,
    Mean_Assigned = mean_assigned,
    Mode_Assigned = mode_assigned,
    mu = mu_mean,
    sigma2 = sigma2_mean,
    pi = pi_mean,
    K = K,
    C = C,
    maxK = maxK,
    category_levels = levels(cat_factor),
    mcmc_samples = list(
      mu = mu_samples,
      sigma2 = sigma2_samples,
      pi = pi_samples,
      z = z_samples
    ),
    data = data,
    varY = varY,
    varCLASS = varCLASS,
    call = match.call()
  )
  
  class(out) = "fuss_COVARIATE_mcmc"
  return(out)
}

# ---- S3 methods ----

#' @export
print.fuss_COVARIATE_mcmc <- function(x, ...) {
  cat("\nUnified Mixture Model with Category-Specific Parameters\n")
  cat("  Categories (C):", x$C, "\n")
  cat("  Components (K):", x$K, "\n")
  cat("\nFirst 6 rows of probability matrix:\n")
  print(head(x$prob_matrix, 6))
  cat("\nFirst 6 assignments:\n")
  print(head(x$Assigned_Group, 6))
  invisible(x)
}

#' @export
summary.fuss_COVARIATE_mcmc <- function(object, ...) {
  cat("\nUnified Mixture Model with Category-Specific Parameters\n")
  cat("  Categories (C):", object$C, "\n")
  cat("  Components (K):", object$K, "\n")
  cat("\nPosterior means (mu):\n")
  print(round(object$mu, 4))
  cat("\nPosterior variances (sigma2):\n")
  print(round(object$sigma2, 4))
  cat("\nPosterior weights (pi):\n")
  print(round(object$pi, 4))
  cat("\nAssignment summary:\n")
  print(table(object$Assigned_Group))
  invisible(object)
}

#' @export
plot.fuss_COVARIATE_mcmc <- function(x, ...) {
  mu_samples = x$mcmc_samples$mu
  C = dim(mu_samples)[1]
  K = dim(mu_samples)[2]
  if (C > 0 && K > 0) {
    old_par = par(no.readonly = TRUE)
    on.exit(par(old_par))
    
    n_rows = min(C, 3)
    n_cols = min(K, 3)
    par(mfrow = c(n_rows, n_cols), mar = c(3, 3, 2, 1))
    
    for (c in 1:n_rows) {
      for (k in 1:n_cols) {
        plot(mu_samples[c, k, ], type = "l",
             main = paste0("mu[", c, ",", k, "]"),
             xlab = "Iteration", ylab = "Value",
             col = "steelblue", lwd = 0.8)
        abline(h = mean(mu_samples[c, k, ]), col = "firebrick2", lwd = 1.5, lty = 2)
      }
    }
  } else {
    message("No MCMC samples to plot.")
  }
}

# ---- MCMC diagnostic functions ----

#' Compute Effective Sample Size (ESS) for MCMC chains
#'
#' @param x Numeric vector or matrix of MCMC samples
#' @return Effective sample size
#' @export
effective_sample_size <- function(x) {
  if (is.vector(x)) {
    n = length(x)
    acf_vals = stats::acf(x, lag.max = min(100, n/2), plot = FALSE)$acf
    rho_sum = sum(acf_vals[-1][acf_vals[-1] > 0])
    ess = n / (1 + 2 * rho_sum)
    return(max(1, ess))
  } else if (is.matrix(x)) {
    return(apply(x, 2, effective_sample_size))
  } else {
    stop("'x' must be a vector or matrix")
  }
}

#' Compute Gelman-Rubin R_hat for MCMC chains
#'
#' @param chains List of either:
#'   - matrices (each row = parameter, each column = iteration), or
#'   - fuss_COVARIATE_mcmc result objects
#' @return Vector of R_hat values (one per parameter)
#' @export
gelman_rubin <- function(chains) {
  if (!is.list(chains) || length(chains) < 2) {
    stop("'chains' must be a list of at least 2 chains")
  }
  
  # ---- Extract matrices from result objects ----
  chain_matrices = list()
  
  for (i in 1:length(chains)) {
    if (inherits(chains[[i]], "fuss_COVARIATE_mcmc")) {
      # Extract mu samples from result object
      cube = chains[[i]]$mcmc_samples$mu
      if (is.null(cube)) {
        stop("Chain ", i, " has no mcmc_samples$mu")
      }
      # Reshape cube to matrix: parameters (C*K) x iterations
      C = dim(cube)[1]
      K = dim(cube)[2]
      n_samples = dim(cube)[3]
      chain_matrices[[i]] = matrix(cube, nrow = C * K, ncol = n_samples)
    } else if (is.matrix(chains[[i]])) {
      # Already a matrix: rows = parameters, cols = iterations
      chain_matrices[[i]] = chains[[i]]
    } else {
      stop("Chain ", i, " must be a fuss_COVARIATE_mcmc object or a matrix")
    }
  }
  
  # ---- Check dimensions ----
  n_params = nrow(chain_matrices[[1]])   # Number of parameters
  n = ncol(chain_matrices[[1]])          # Number of iterations per chain
  
  for (i in 1:length(chain_matrices)) {
    if (nrow(chain_matrices[[i]]) != n_params || ncol(chain_matrices[[i]]) != n) {
      stop("All chains must have the same dimensions")
    }
  }
  
  n_chains = length(chain_matrices)
  
  # ---- Calculate chain means ----
  # Each column is a chain, each row is a parameter
  chain_means = matrix(NA, nrow = n_params, ncol = n_chains)
  for (j in 1:n_chains) {
    chain_means[, j] = rowMeans(chain_matrices[[j]])
  }
  
  # Overall mean across all chains (per parameter)
  overall_mean = rowMeans(chain_means)
  
  # ---- Between-chain variance (B) ----
  # B = n/(m-1) * sum_j (mean_j - mean_overall)^2
  B = rep(0, n_params)
  for (j in 1:n_chains) {
    diff = chain_means[, j] - overall_mean
    B = B + diff^2
  }
  B = B * n / (n_chains - 1)
  
  # ---- Within-chain variance (W) ----
  # W = 1/m * sum_j var(chain_j)
  W = rep(0, n_params)
  for (j in 1:n_chains) {
    W = W + apply(chain_matrices[[j]], 1, var)
  }
  W = W / n_chains
  
  # ---- R_hat ----
  # V_hat = (n-1)/n * W + 1/n * B
  V_hat = (n - 1) / n * W + B / n
  R_hat = sqrt(V_hat / W)
  
  names(R_hat) = paste0("param_", 1:n_params)
  return(R_hat)
}

#' Check MCMC convergence
#'
#' @param object A \code{fuss_COVARIATE_mcmc} object
#' @param chains Optional list of multiple chains for Gelman-Rubin
#' @return Convergence diagnostics
#' @export
check_convergence <- function(object, chains = NULL) {
  if (!inherits(object, "fuss_COVARIATE_mcmc")) {
    stop("'object' must be of class 'fuss_COVARIATE_mcmc'")
  }
  
  if (is.null(chains)) {
    # Single chain: ESS
    mu_samples = object$mcmc_samples$mu
    C = dim(mu_samples)[1]
    K = dim(mu_samples)[2]
    
    ess = matrix(NA, C, K)
    for (c in 1:C) {
      for (k in 1:K) {
        ess[c, k] = effective_sample_size(mu_samples[c, k, ])
      }
    }
    rownames(ess) = paste0("Cat_", 1:C)
    colnames(ess) = paste0("Comp_", 1:K)
    
    cat("\nEffective Sample Size (ESS):\n")
    print(round(ess, 0))
    cat("\nESS should be > 400 for reliable inference.\n")
    return(invisible(ess))
    
  } else {
    # Multiple chains: Gelman-Rubin
    chain_matrices = lapply(chains, function(ch) {
      if (!inherits(ch, "fuss_COVARIATE_mcmc")) {
        stop("Each element of 'chains' must be a 'fuss_COVARIATE_mcmc' object")
      }
      cube = ch$mcmc_samples$mu
      C = dim(cube)[1]
      K = dim(cube)[2]
      matrix(cube, nrow = C * K, ncol = dim(cube)[3])
    })
    
    R_hat = gelman_rubin(chain_matrices)
    names(R_hat) = paste0("mu[", rep(1:C, each = K), ",", rep(1:K, C), "]")
    
    cat("\nGelman-Rubin R_hat:\n")
    print(round(R_hat, 3))
    cat("\nR_hat should be < 1.1 for convergence.\n")
    return(invisible(R_hat))
  }
}