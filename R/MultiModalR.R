# MULTIMODALR - Fast Bayesian Probability Estimation for Multimodal Categorical Data
# Version: 1.0.0
# Speed-optimized MCMC implementation (Metropolis-Hastings-within-partial-Gibbs)
# Based on MINLAM (depreciated) by Gergő Diószegi

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

#' Obtain the number of groups based on bandwidth selection with SJ tuning
#'
#' @param y vector, input data of a distribution
#' @param sj_adjust adjustment factor for sj-dpi bandwidth (default: 1.0)
#' @return data frame with methods and number of groups
#' @export
get_NGRP <- function(y, sj_adjust = 1.0) {
  
  # Calculate bandwidths
  bwRT = stats::bw.nrd(y)                 # Scott's 1.06
  bwBCV = stats::bw.bcv(y)                # Biased cross-validation
  
  # Sheather & Jones' direct plug-in method with tuning
  bwSJ_base = stats::bw.SJ(y, method = "dpi")
  bwSJ_tuned = bwSJ_base * sj_adjust
  
  # Ensure reasonable bounds
  data_range = diff(range(y))
  min_bw = data_range / 100  # Minimum: 1% of range
  max_bw = data_range / 3    # Maximum: 33% of range
  bwSJ_tuned = max(min(bwSJ_tuned, max_bw), min_bw)
  
  # Detect modes
  nrd_nmod = multimode::nmodes(y, bw = bwRT)
  bcv_nmod = multimode::nmodes(y, bw = bwBCV)
  sj_nmod = multimode::nmodes(y, bw = bwSJ_tuned)
  
  return(data.frame(
    Method = c("nrd", "bcv", "sj-dpi"),
    Bandwidth = c(bwRT, bwBCV, bwSJ_tuned),
    n_grp = c(nrd_nmod, bcv_nmod, sj_nmod)
  ))
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
#' @param within numeric, range for grouping modes (default: 0.1)
#' @return data frame with grouped modes
#' @export
group_MODES <- function(df, within = 0.1) {
  df = df %>%
    dplyr::arrange(Est_Mode) %>%
    dplyr::mutate(group = cumsum(c(1, diff(Est_Mode) > within)))
  
  df %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(Est_Mode = mean(Est_Mode), .groups = "drop")
}

#' Fast MCMC for mixture models (C++ implementation)
#' 
#' @param y Numeric vector of data
#' @param grp Number of mixture components
#' @param prior_means Prior means for components (optional)
#' @param ids Vector of IDs for validation (required)
#' @param n_iter Number of MCMC iterations (default: 1000)
#' @param burnin Burn-in period (default: 500)
#' @param proposal_sd Proposal standard deviation for component means (default: 0.15)
#' @param seed Random seed
#' @return List with MCMC results
#' @export
MM_MH <- function(y, grp, prior_means = NULL, ids,
                  n_iter = 1000, burnin = 500,
                  proposal_sd = 0.15, seed = NULL) {
  
  # Validate IDs - they are required
  if(missing(ids)) {
    stop("'ids' argument is required for validation.")
  }
  
  if(is.null(ids)) {
    stop("'ids' cannot be NULL. Provide a vector of IDs for validation.")
  }
  
  if(length(ids) != length(y)) {
    stop("ID length (", length(ids), ") doesn't match y length (", length(y), ").")
  }
  
  if(any(is.na(ids))) {
    warning("NA values found in IDs. These will be preserved but may cause issues in validation.")
  }
  
  if(any(duplicated(ids))) {
    warning("Duplicate IDs found. This may cause issues in validation.")
  }
  
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
  
  # Store the IDs in the result
  result$ids = ids
  
  return(result)
}

#' Create output data frame
#' 
#' Converts MCMC results to exact CSV format
#' 
#' @param mcmc_result Output from MM_MH() - MUST contain ids in mcmc_result$ids if needed
#' @param y_original Original y values (if different from mcmc_result$y)
#' @param group_original Original group assignments (optional)
#' @param main_class Category/class name
#' @param max_groups Maximum number of groups for output columns
#' @return Data frame in MINLAM CSV format
#' @export
create_MM_output <- function(mcmc_result, y_original = NULL, 
                             group_original = NULL, 
                             main_class = "", max_groups = 5) {
  
  n = length(mcmc_result$y)
  n_components = mcmc_result$n_components
  
  # Get IDs from mcmc_result (should have been added by MM_MH if provided)
  ids = mcmc_result$ids
  
  # Check ID length if IDs exist
  if(!is.null(ids) && length(ids) != n) {
    warning("ID length in mcmc_result (", length(ids), ") doesn't match data length (", n, "). IDs will be omitted.")
    ids = NULL
  }
  
  # Use original y if provided
  if(!is.null(y_original)) {
    if(length(y_original) != n) {
      warning("y_original length doesn't match. Using MCMC y values.")
      y_values = mcmc_result$y
    } else {
      y_values = y_original
    }
  } else {
    y_values = mcmc_result$y
  }
  
  # Create the data frame with or without ID
  if(!is.null(ids)) {
    df = data.frame(
      ID = ids,
      y = y_values,
      stringsAsFactors = FALSE
    )
  } else {
    df = data.frame(
      y = y_values,
      stringsAsFactors = FALSE
    )
  }
  
  # Add original group if provided
  if(!is.null(group_original)) {
    if(length(group_original) == n) {
      df$Group = group_original
    } else {
      warning("Length of group_original doesn't match data. Skipping.")
    }
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

#' Core function for mixture model MCMC with ID support
#' 
#' @param data Input data frame
#' @param varCLASS Character, category variable name (required)
#' @param varY Character, value variable name (required)
#' @param varID Character, ID variable name (required)
#' @param method Density estimator method ("sj-dpi", "bcv", "nrd")
#' @param within Range parameter for mode grouping (default: 0.1)
#' @param maxNGROUP Maximum number of groups
#' @param out_dir Output directory for CSV files (if NULL, returns data frame)
#' @param n_iter Number of MCMC iterations (default: 1000)
#' @param burnin Burn-in period (default: 500)
#' @param proposal_sd Proposal standard deviation for component means (default: 0.15)
#' @param sj_adjust Adjustment factor for sj-dpi bandwidth detector (default: 1.0)
#' @return Data frame or writes CSV files to out_dir
#' @export
get_PROBCLASS_MH <- function(data, varCLASS, varY, varID, 
                             method = "sj-dpi", within = 0.1, maxNGROUP = 5, 
                             out_dir = NULL, n_iter = 1000,
                             burnin = 500, proposal_sd = 0.15,
                             sj_adjust = 1.0) {
  
  # Validate inputs
  if(!is.data.frame(data)) {
    stop("data must be a data frame")
  }
  
  if(!varCLASS %in% names(data)) {
    stop("varCLASS '", varCLASS, "' not found in data")
  }
  
  if(!varY %in% names(data)) {
    stop("varY '", varY, "' not found in data")
  }
  
  if(!varID %in% names(data)) {
    stop("varID '", varID, "' not found in data")
  }
  
  # Warn if sj_adjust is used with non-sj-dpi method
  if(method != "sj-dpi" && sj_adjust != 1.0) {
    warning("sj_adjust = ", sj_adjust, " is IGNORED for method = '", method, "'.\n",
            "  sj_adjust only affects 'sj-dpi' (Sheather-Jones) method.\n",
            "  Using sj_adjust = 1.0 for method = '", method, "'.")
  }
  
  categories = unique(data[[varCLASS]])
  results = list()
  
  for(cat in categories) {
    cat_data = data[data[[varCLASS]] == cat, ]
    y = cat_data[[varY]]
    
    # Extract IDs
    ids = cat_data[[varID]]
    
    message("Processing category: ", cat)
    
    # Determine effective sj_adjust value
    # Only use sj_adjust for method = "sj-dpi", otherwise use 1.0
    effective_sj_adjust <- ifelse(method == "sj-dpi", sj_adjust, 1.0)
    
    if(method == "sj-dpi") {
      message("  Using SJ-dpi with adjustment: ", sj_adjust)
    } else {
      message("  Using method: ", method, " (sj_adjust ignored)")
    }
    
    # Mode detection WITH SJ TUNING
    n_grp_df = get_NGRP(y, sj_adjust = effective_sj_adjust)
    n_grp = n_grp_df %>% 
      dplyr::filter(Method == method) %>% 
      dplyr::pull(n_grp)
    n_grp = min(max(n_grp, 2), maxNGROUP)
    
    # Get bandwidth used for debugging
    bandwidth_used = n_grp_df %>% 
      dplyr::filter(Method == method) %>% 
      dplyr::pull(Bandwidth)
    
    message("  Detected ", n_grp, " components (BW = ", round(bandwidth_used, 4), ")")
    
    # Get mode locations
    modes_df = get_MODES(y, nmod = n_grp)
    modes_grouped = group_MODES(modes_df, within = within)
    
    n_components = nrow(modes_grouped)
    prior_means = modes_grouped$Est_Mode
    
    message("  After grouping: ", n_components, " components")
    
    # Run MCMC WITH IDs
    mh_result = MM_MH(
      y = y,
      grp = n_components,
      prior_means = prior_means,
      ids = ids,  
      n_iter = n_iter,
      burnin = burnin,
      proposal_sd = proposal_sd,
      seed = 123
    )
    
    message("  Mean acceptance: ", paste(round(mh_result$mean_acceptance, 3), collapse = ", "))
    
    # Create output - IDs are already in mcmc_result from MM_MH
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
      message("  Written: ", filepath)
    }
  }
  
  # Combine and return results if no out_dir
  if(is.null(out_dir)) {
    combined_result = do.call(rbind, results)
    message("Analysis complete. Result has ", nrow(combined_result), " rows.")
    
    return(combined_result)
  } else {
    message("Results written to: ", out_dir)
    return(invisible(NULL))
  }
}

#' Parallel wrapper around get_PROBCLASS_MH
#' 
#' @param data Input data frame
#' @param varCLASS Character, category variable name (required)
#' @param varY Character, value variable name (required)
#' @param varID Character, ID variable name (required)
#' @param method Density estimator method ("sj-dpi", "bcv", "nrd") (default: "sj-dpi")
#' @param within Range parameter
#' @param maxNGROUP Maximum number of groups
#' @param out_dir Output directory for CSV files (if NULL, returns combined data frame)
#' @param n_workers Number of parallel workers
#' @param n_iter Number of MCMC iterations (default: 1000)
#' @param burnin Burn-in period (default: 500)
#' @param proposal_sd Proposal standard deviation for component means (default: 0.15)
#' @param sj_adjust Adjustment factor for sj-dpi bandwidth (default: 1.0, smaller -> more modes, higher -> fewer modes)
#' @return Data frame (if out_dir is NULL) or writes CSV files
#' @export
fuss_PARALLEL <- function(data, varCLASS, varY, varID, 
                          method = "sj-dpi", within = 1, maxNGROUP = 5, 
                          out_dir = NULL, n_workers = 4,  
                          n_iter = 1000, burnin = 500,
                          proposal_sd = 0.15, sj_adjust = 1.0) {
  
  # Validate inputs
  if(!is.data.frame(data)) {
    stop("data must be a data frame")
  }
  
  if(!varCLASS %in% names(data)) {
    stop("varCLASS '", varCLASS, "' not found in data")
  }
  
  if(!varY %in% names(data)) {
    stop("varY '", varY, "' not found in data")
  }
  
  if(!varID %in% names(data)) {
    stop("varID '", varID, "' not found in data")
  }
  
  # Split data by category
  categories = unique(data[[varCLASS]])
  data_list = split(data, data[[varCLASS]])
  
  message("Processing ", length(categories), " categories in parallel with ", 
          n_workers, " workers")
  
  # Warn if sj_adjust is used with non-sj-dpi method
  if(method != "sj-dpi" && sj_adjust != 1.0) {
    warning("sj_adjust = ", sj_adjust, " is IGNORED for method = '", method, "'.\n",
            "  sj_adjust only affects 'sj-dpi' (Sheather-Jones) method.\n",
            "  Using sj_adjust = 1.0 for method = '", method, "'.")
  }
  
  # Setup parallel processing
  future::plan(future::multisession, workers = n_workers)
  
  # Process each category in parallel
  result_list = furrr::future_map(data_list, function(cat_data) {
    
    # Get category name
    cat_name = as.character(unique(cat_data[[varCLASS]])[1])
    message("Processing category: ", cat_name)
    
    # Determine effective sj_adjust value
    effective_sj_adjust <- ifelse(method == "sj-dpi", sj_adjust, 1.0)
    
    # Use get_PROBCLASS_MH with sj_adjust parameter
    result = get_PROBCLASS_MH(
      data = cat_data,
      varCLASS = varCLASS,
      varY = varY,
      varID = varID,
      method = method,
      within = within,
      maxNGROUP = maxNGROUP,
      out_dir = out_dir,
      n_iter = n_iter,
      burnin = burnin,
      proposal_sd = proposal_sd,
      sj_adjust = effective_sj_adjust
    )
    
    return(result)
    
  }, .options = furrr::furrr_options(
    packages = "MultiModalR",
    seed = TRUE
  ), .progress = TRUE)
  
  # Reset to sequential
  future::plan(future::sequential)
  
  # If out_dir is NULL, combine and return results
  if(is.null(out_dir)) {
    # Remove NULL results (happens when out_dir is provided in get_PROBCLASS_MH)
    result_list = result_list[!sapply(result_list, is.null)]
    
    if(length(result_list) > 0) {
      combined_result = do.call(rbind, result_list)
      message("Parallel analysis complete. Combined result has ", 
              nrow(combined_result), " rows.")
      
      return(combined_result)
    } else {
      message("All results written to output directory.")
      return(invisible(NULL))
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
#' @param id_col Character, name of the ID column in observed_df (default: "ID")
#' @param pattern Pattern to match CSV files (default: "^df_")
#' @return ggplot object showing validation results
#' @export
plot_VALIDATION <- function(csv_dir, observed_df, 
                            subpop_col = "Subpopulation", 
                            value_col = "Value",
                            id_col = "ID",
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
  predicted_list = list()
  for(f in FIL) {
    predicted_list[[f]] = readr::read_csv(f, show_col_types = FALSE) %>%
      as.data.frame() %>%
      dplyr::mutate(Main_Class = as.character(Main_Class))
  }
  
  predicted = purrr::map_dfr(predicted_list, ~ .x) %>%
    dplyr::mutate(Main_Class = gsub("__", "::", Main_Class))
  
  # Extract true subgroup numbers
  observed_df[[subpop_col]] = as.numeric(
    gsub("Group ", "", as.character(observed_df[[subpop_col]]))
  )
  
  # Ensure ID columns are character type
  observed_df[[id_col]] = as.character(observed_df[[id_col]])
  predicted[[id_col]] = as.character(predicted[[id_col]])
  
  # Match by ID only
  matched = dplyr::inner_join(observed_df, predicted, by = id_col)
  
  # Calculate overall accuracy
  matched$correct = matched[[subpop_col]] == matched$Assigned_Group
  overall_accuracy = round(mean(matched$correct, na.rm = TRUE) * 100, 1)
  
  # Calculate accuracy per Main_Class
  accuracy_by_class = matched %>%
    group_by(Main_Class) %>%
    summarize(
      accuracy = round(mean(correct, na.rm = TRUE) * 100, 1),
      n = n(),
      .groups = 'drop'
    )
  
  # Create label data for plot
  label_data = accuracy_by_class %>%
    mutate(Percent = format(accuracy, nsmall = 1))
  
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