# INLAMM - Fast Bayesian Probability Estimation for Multimodal Categorical Data
# Version: 1.0.0
# Speed-optimized MCMC implementation (Metropolis-Hastings-within-partial-Gibbs)
# Based on MINLAM by DijoG

#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib MultiModalR, .registration = TRUE
NULL

#' Check and install required packages
#' 
#' @return Installs missing packages
#' @export
check_PACKS <- function() {
  required_packages = c(
    "tidyverse", "furrr", "future", "multimode", "Rcpp", "RcppArmadillo",
    "truncnorm", "tictoc", "dplyr", "tidyr", "purrr", "readr", "ggplot2"
  )
  
  missing_packages = required_packages[!required_packages %in% installed.packages()]
  
  if(length(missing_packages) > 0) {
    message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
    install.packages(missing_packages)
  } else {
    message("All required packages are installed.")
  }
}

#' Obtain the number of groups based on bandwidth selection
#'
#' @param y vector, input data of a distribution
#' @return data frame with methods and number of groups
#' @export
get_NGRP <- function(y) {
  bwRT = stats::bw.nrd(y)                 # Scott's 1.06
  bwSJ = stats::bw.SJ(y, method = "dpi")  # Sheather & Jones' direct plug-in method
  bwBCV = stats::bw.bcv(y)                # Biased cross-validation
  
  nrd_nmod = multimode::nmodes(y, bw = bwRT)
  sj_nmod = multimode::nmodes(y, bw = bwSJ)
  bcv_nmod = multimode::nmodes(y, bw = bwBCV)
  
  return(data.frame(Method = c("nrd", "bcv", "dpi"),
                    n_grp = c(nrd_nmod, bcv_nmod, sj_nmod)))
}

#' Extract mode statistics from a mode forest
#'
#' @param y vector, input data of a distribution
#' @param nmod numeric, count/number of subpopulations/subgroups  
#' @return data frame with estimated modes
#' @export
get_MODES <- function(y, nmod) {
  mm = multimode::locmodes(y, mod0 = nmod)
  modes = mm$locations[seq(1, length(mm$locations), by = 2)]
  mode_df = data.frame(Est_Mode = modes,
                       Group = 1:length(modes))
  return(mode_df)
}

#' Group/merge modes if they are within a given range
#'
#' @param df data frame containing samples of a distribution in its 'Est_Mode' variable
#' @param within numeric, range for grouping modes
#' @return data frame with grouped modes
#' @export
group_MODES <- function(df, within = 0.03) {
  df = df %>%
    dplyr::arrange(Est_Mode) %>%
    dplyr::mutate(group = cumsum(c(1, diff(Est_Mode) > within)))
  
  df %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(Est_Mode = mean(Est_Mode), .groups = "drop")
}

#' Fast MCMC for mixture models (C++ implementation)
#' 
#' Fixed version that prevents empty groups
#' 
#' @param y Numeric vector of data
#' @param grp Number of mixture components
#' @param prior_means Prior means for components (optional)
#' @param n_iter Number of MCMC iterations (default: 1000)
#' @param burnin Burn-in period (default: 500)
#' @param proposal_sd Proposal standard deviation for component means (default: 0.15)
#' @param seed Random seed
#' @return List with MCMC results
#' @export
MM_MH <- function(y, grp, prior_means = NULL,
                  n_iter = 1000, burnin = 500,
                  proposal_sd = 0.15, seed = NULL) {
  
  if(is.null(seed)) {
    seed = sample.int(.Machine$integer.max, 1)
  }
  
  if(is.null(prior_means)) {
    prior_means = seq(min(y), max(y), length.out = grp)
  }
  
  # Ensure prior_means matches grp
  if(length(prior_means) != grp) {
    warning("prior_means length doesn't match grp. Using equally spaced means.")
    prior_means = seq(min(y), max(y), length.out = grp)
  }
  
  result = MM_MH_cpp(
    y = as.numeric(y),
    prior_means = as.numeric(prior_means),
    n_iter = as.integer(n_iter),
    burnin = as.integer(burnin),
    proposal_sd = as.numeric(proposal_sd),
    seed = as.integer(seed)
  )
  
  return(result)
}

#' Create output data frame
#' 
#' Converts MCMC results to exact CSV format
#' 
#' @param mcmc_result Output from MM_MH()
#' @param y_original Original y values (if different)
#' @param group_original Original group assignments (optional)
#' @param main_class Category/class name
#' @param max_groups Maximum number of groups for output columns
#' @return Data frame in MINLAM CSV format
#' @export
create_MM_output <- function(mcmc_result, y_original = NULL, 
                             group_original = NULL, main_class = "",
                             max_groups = 5) {
  
  n = length(mcmc_result$y)
  n_components = mcmc_result$n_components
  
  # Use original y if provided
  if(!is.null(y_original)) {
    y_values = y_original
  } else {
    y_values = mcmc_result$y
  }
  
  # Create the data frame
  df = data.frame(
    y = y_values,
    stringsAsFactors = FALSE
  )
  
  # Add original group if provided
  if(!is.null(group_original)) {
    df$Group = group_original
  }
  
  # Add probability columns for actual components
  for(k in 1:n_components) {
    df[[paste0("Group_", k)]] = mcmc_result$z_probs[, k]
  }
  
  # Add NA columns for missing groups (up to max_groups)
  for(k in (n_components + 1):max_groups) {
    df[[paste0("Group_", k)]] = NA_real_
  }
  
  # Add assignment columns
  df$Assigned_Group = mcmc_result$assigned_group
  df$Min_Assigned = mcmc_result$min_assigned
  df$Max_Assigned = mcmc_result$max_assigned
  df$Mean_Assigned = mcmc_result$mean_assigned
  df$Mode_Assigned = mcmc_result$mode_assigned
  
  # Ensure Main_Class is character (not factor or logical)
  df$Main_Class = as.character(main_class)
  
  return(df)
}

#' Fast version of get_PROBCLASS_MH
#' 
#' @param data Input data frame
#' @param varCLASS Character, category variable name
#' @param varY Character, value variable name
#' @param method Density estimator method (for compatibility)
#' @param within Range parameter (for compatibility)
#' @param maxNGROUP Maximum number of groups
#' @param out_dir Output directory for CSV files (if NULL, returns data frame)
#' @param n_iter Number of MCMC iterations (default: 1000)
#' @param burnin Burn-in period (default: 500)
#' @param proposal_sd Proposal standard deviation for component means (default: 0.15)
#' @return Data frame or writes CSV files to out_dir
#' @export
get_PROBCLASS_MH <- function(data, varCLASS, varY, method = "dpi", 
                             within = 0.03, maxNGROUP = 5, 
                             out_dir = NULL,  n_iter = 1000,
                             burnin = 500, proposal_sd = 0.15) {
  
  categories = unique(data[[varCLASS]])
  results = list()
  
  for(cat in categories) {
    cat_data = data[data[[varCLASS]] == cat, ]
    y = cat_data[[varY]]
    
    # Mode detection
    n_grp_df = get_NGRP(y)
    n_grp = n_grp_df %>% 
      dplyr::filter(Method == method) %>% 
      dplyr::pull(n_grp)
    n_grp = min(max(n_grp, 3), maxNGROUP)
    
    message("Category ", cat, ": method ", method, " detected ", n_grp, " components")
    
    # Get mode locations
    modes_df = get_MODES(y, nmod = n_grp)
    modes_grouped = group_MODES(modes_df, within = within)
    
    n_components = nrow(modes_grouped)
    prior_means = modes_grouped$Est_Mode
    
    message("  After grouping: ", n_components, " components with means: ", 
            paste(round(prior_means, 2), collapse = ", "))
    
    # Run MCMC
    mh_result = MM_MH(
      y = y,
      grp = n_components,
      prior_means = prior_means,
      n_iter = 1000,
      burnin = 500,
      proposal_sd = 0.15,
      seed = 123
    )
    
    message("  Mean acceptance: ", paste(round(mh_result$mean_acceptance, 3), collapse = ", "))
    
    # Create output
    output_df = create_MM_output(
      mcmc_result = mh_result,
      y_original = y,
      group_original = NULL,
      main_class = as.character(cat),
      max_groups = maxNGROUP
    )
    
    results[[as.character(cat)]] = output_df
    
    # Write to CSV if out_dir provided
    if(!is.null(out_dir)) {
      if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
      filename = paste0("df_", cat, ".csv")
      filepath = file.path(out_dir, filename)
      write.csv(output_df, filepath, row.names = FALSE, quote = TRUE)
      message("Written: ", filepath)
    }
  }
  
  # Combine and return results if no out_dir
  if(is.null(out_dir)) {
    return(do.call(rbind, results))
  } else {
    message("Results written to: ", out_dir)
    return(invisible(NULL))
  }
}

#' Fast parallel version
#' 
#' @param data List of data frames (grouped by GROUP variable)
#' @param varCLASS Character, category variable name
#' @param varY Character, value variable name
#' @param method Density estimator method
#' @param within Range parameter
#' @param maxNGROUP Maximum number of groups
#' @param out_dir Output directory for CSV files (if NULL, returns combined data frame)
#' @param n_workers Number of parallel workers
#' @param n_iter Number of MCMC iterations (default: 1000)
#' @param burnin Burn-in period (default: 500)
#' @param proposal_sd Proposal standard deviation for component means (default: 0.15)
#' @return Data frame (if out_dir is NULL) or writes CSV files
#' @export
fuss_PARALLEL <- function(data, varCLASS, varY, method = "dpi", 
                          within = 1, maxNGROUP = 5, 
                          out_dir = NULL, n_workers = 4,  
                          n_iter = 1000, burnin = 500,
                          proposal_sd = 0.15) {
  
  future::plan(future::multisession, workers = n_workers)
  
  result_list = furrr::future_map(data, function(.x) {
    get_PROBCLASS_MH(.x, varCLASS, varY, method,
                     within, maxNGROUP, out_dir)
  }, .options = furrr::furrr_options(
    packages = "INLAMM",
    seed = TRUE
  ), .progress = TRUE)
  
  # If out_dir is NULL, combine and return results
  if(is.null(out_dir)) {
    # Remove NULL results (from files that were written to disk)
    result_list = result_list[!sapply(result_list, is.null)]
    if(length(result_list) > 0) {
      return(dplyr::bind_rows(result_list))
    } else {
      return(NULL)
    }
  } else {
    message("All results written to: ", out_dir)
    return(invisible(NULL))
  }
}

#' Plot validation of subgroup assignments (handles both balanced and imbalanced data)
#'
#' @param csv_dir Directory containing CSV files from create_MM_output
#' @param observed_df Original data frame with true subgroups
#' @param subpop_col Character, name of the true subgroup column in observed_df (default: "Subpopulation")
#' @param value_col Character, name of the value column in observed_df (default: "Value")
#' @param pattern Pattern to match CSV files (default: "^df_")
#' @return ggplot object showing validation results
#' @export
plot_VALIDATION <- function(csv_dir, observed_df, 
                            subpop_col = "Subpopulation", 
                            value_col = "Value",
                            pattern = "^df_") {
  
  # Check required packages
  required = c("ggplot2", "dplyr", "readr", "purrr", "tidyr")
  missing = required[!required %in% installed.packages()]
  if(length(missing) > 0) {
    stop("Missing packages: ", paste(missing, collapse = ", "),
         "\nPlease install with: install.packages(c('", paste(missing, collapse = "', '"), "'))")
  }
  
  # Get list of CSV files
  FIL = list.files(csv_dir, pattern = pattern, full.names = TRUE) 
  
  if(length(FIL) == 0) {
    stop("No CSV files found in ", csv_dir, " matching pattern: ", pattern)
  }
  
  # Read and combine predicted data
  predicted = purrr::map_dfr(FIL, ~ readr::read_csv(.x, 
                                                    col_types = "dddddddddddc", 
                                                    show_col_types = FALSE)) %>%
    as.data.frame() %>%
    dplyr::mutate(Main_Class = factor(as.character(Main_Class))) %>%
    dplyr::arrange(y)
  
  # Extract true subgroup numbers (remove "Group " prefix if present)
  observed_df[[subpop_col]] = as.numeric(
    gsub("Group ", "", as.character(observed_df[[subpop_col]]))
  )
  
  # Sort observed data by value
  observed_df = observed_df %>% 
    dplyr::arrange(.data[[value_col]])
  
  # Debug information
  message("Dataset sizes:")
  message("  Predicted: ", nrow(predicted), " rows")
  message("  Observed:  ", nrow(observed_df), " rows")
  
  # Check alignment
  if(nrow(observed_df) != nrow(predicted)) {
    warning("Observed (n=", nrow(observed_df), ") and predicted (n=", nrow(predicted), 
            ") have different numbers of observations. Aligning by sorted values.")
  }
  
  # Take the smaller number of observations to compare
  n_compare = min(nrow(predicted), nrow(observed_df))
  
  if(n_compare == 0) {
    stop("No observations to compare.")
  }
  
  # Compare first n_compare observations (assumes sorted by value)
  matching_indices = which(
    observed_df[[subpop_col]][1:n_compare] == predicted$Assigned_Group[1:n_compare]
  )
  
  # Calculate overall accuracy
  overall_accuracy = round(length(matching_indices) / n_compare * 100, 1)
  message("Overall matching accuracy: ", overall_accuracy, "% (based on ", n_compare, " aligned observations)")
  
  # Calculate accuracy per Main_Class (handles imbalanced data automatically)
  main_classes = unique(predicted$Main_Class)
  accuracy_list = list()
  
  for(main_class in main_classes) {
    pred_idx = which(predicted$Main_Class == main_class)
    
    if(length(pred_idx) == 0) {
      accuracy_list[[as.character(main_class)]] = 0
      next
    }
    
    # Find corresponding observed indices (assumes same order after sorting)
    if(length(pred_idx) <= n_compare) {
      matches = sum(
        observed_df[[subpop_col]][pred_idx] == predicted$Assigned_Group[pred_idx],
        na.rm = TRUE
      )
      accuracy = ifelse(length(pred_idx) > 0, 
                        round(matches / length(pred_idx) * 100, 1), 
                        0)
    } else {
      accuracy = 0
    }
    
    accuracy_list[[as.character(main_class)]] = accuracy
    message("  ", main_class, ": ", accuracy, "%")
  }
  
  # Create label data for plot
  if(length(accuracy_list) > 0) {
    label_data = data.frame(
      Main_Class = names(accuracy_list),
      Percent = format(round(unlist(accuracy_list), 1), nsmall = 1),
      stringsAsFactors = FALSE
    )
  } else {
    label_data = data.frame(
      Main_Class = as.character(main_classes),
      Percent = "0.0",
      stringsAsFactors = FALSE
    )
  }
  
  # Create color palette for subgroups
  n_groups = length(unique(predicted$Assigned_Group))
  group_colors = c("firebrick2", "forestgreen", "cyan3", "gold", "purple", "orange")
  group_colors = group_colors[1:min(n_groups, length(group_colors))]
  
  # Plot validation
  p = ggplot2::ggplot(predicted, ggplot2::aes(x = y)) +
    ggplot2::geom_density(col = NA, fill = "grey98", adjust = 0.8) +
    ggplot2::geom_jitter(ggplot2::aes(y = 0.05, color = factor(Assigned_Group)), 
                         height = 0.05, size = 2, shape = 16, alpha = .5) + 
    ggplot2::scale_color_manual(values = group_colors, 
                                name = "Assigned Groups") +  
    ggplot2::theme_dark() +
    ggplot2::labs(title = paste0("Validation of Subgroup Assignments (", overall_accuracy, "% overall)"),
                  x = "Value", y = "Density") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::facet_wrap(~ Main_Class, ncol = 3) +  
    ggplot2::geom_text(data = label_data, ggplot2::aes(x = Inf, y = Inf, 
                                                       label = paste0(Percent, "%")), 
                       hjust = 1.2, vjust = 1.2, size = 5, fontface = "bold", 
                       inherit.aes = FALSE, col = "grey15") +  
    ggplot2::theme(legend.position = "top",
                   legend.key = ggplot2::element_rect(fill = "transparent", color = NA),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(hjust = .5),
                   plot.subtitle = ggplot2::element_text(hjust = .5, size = 10)) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = .7)))
  
  return(p)
}
