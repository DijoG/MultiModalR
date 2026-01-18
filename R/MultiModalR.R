# MULTIMODALR - Fast Bayesian Probability Estimation for Multimodal Categorical Data
# Version: 1.0.0
# Speed-optimized MCMC implementation (Metropolis-Hastings-within-partial-Gibbs and Dirichlet-Multinomial)
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

#' Density height-aware mode detection
#' 
#' Returns mode estimates from FOUR different bandwidth methods.
#' Each method may detect different numbers and locations of modes.
#' 
#' @param y Numeric vector
#' @param adjust Bandwidth adjustment factor (affects "SJ", "nrd", "bcv" methods)
#' @param threshold Relative threshold for significant peaks
#' @return List with mode estimates from multiple methods including density heights
#' @export
get_MODES_enhanced <- function(y, adjust = 1, threshold = 1.0) {
  
  # Estimate densities with different bandwidths
  dSJ = density(y, bw = "SJ", adjust = adjust)
  dNRD = density(y, bw = "nrd", adjust = adjust)
  dBCV = density(y, bw = "bcv", adjust = adjust)
  dUCV = density(y, bw = "ucv")
  
  # Function to find peaks with heights
  find_peaks_with_heights <- function(dens, method_name) {
    
    # Find local maxima
    y_diff1 = diff(dens$y)
    y_diff2 = diff(sign(y_diff1))
    peak_idx = which(y_diff2 < 0) + 1
    
    peaks = dens$x[peak_idx]
    heights = dens$y[peak_idx]
    
    # Filter peaks above threshold
    significant = heights > mean(dens$y) * threshold
    sig_peaks = peaks[significant]
    sig_heights = heights[significant]
    
    # Return as data frame with heights
    if(length(sig_peaks) > 0) {
      data.frame(
        Est_Mode = sig_peaks,
        Density_Height = sig_heights,
        Group = 1:length(sig_peaks),
        Method = method_name
      )
    } else {
      data.frame(
        Est_Mode = numeric(0),
        Density_Height = numeric(0),
        Group = integer(0),
        Method = character(0)
      )
    }
  }
  
  # Get peaks for each method
  peaks_SJ = find_peaks_with_heights(dSJ, "sj-dpi")
  peaks_NRD = find_peaks_with_heights(dNRD, "nrd")
  peaks_BCV = find_peaks_with_heights(dBCV, "bcv")
  peaks_UCV = find_peaks_with_heights(dUCV, "ucv")
  
  # Store bandwidths
  bandwidths = c(
    "sj-dpi" = dSJ$bw,
    "nrd" = dNRD$bw,
    "bcv" = dBCV$bw,
    "ucv" = dUCV$bw
  )
  
  return(list(
    "sj-dpi" = peaks_SJ,
    "nrd" = peaks_NRD,
    "bcv" = peaks_BCV,
    "ucv" = peaks_UCV,
    bandwidths = bandwidths,
    threshold_used = threshold,
    adjust_used = adjust,
    n_obs = length(y)
  ))
}

#' Density height-aware mode grouping
#'
#' @param df data frame containing 'Est_Mode' and 'Density_Height' columns
#' @param within numeric, range for grouping modes (default: 0.1)
#' @return data frame with grouped modes
#' @export
group_MODES_enhanced <- function(df, within = 0.1) {
  
  # Check if we have data
  if(nrow(df) == 0) {
    return(data.frame(
      group = integer(0),
      Est_Mode = numeric(0)
    ))
  }
  
  # Check if density heights are available
  if("Density_Height" %in% names(df)) {
    
    # Enhanced grouping: sort by height (descending)
    modes_df = df %>%
      dplyr::arrange(dplyr::desc(Density_Height)) %>%
      dplyr::mutate(
        processed = FALSE,
        group = NA_integer_
      )
    
    group_counter = 1
    n_modes = nrow(modes_df)
    
    # Process modes from highest to lowest density
    for(i in 1:n_modes) {
      if(!modes_df$processed[i]) {
        current_mode = modes_df$Est_Mode[i]
        
        # Start a new group with this mode
        modes_df$group[i] = group_counter
        modes_df$processed[i] = TRUE
        
        # Find all unprocessed modes within range
        if(i < n_modes) {
          for(j in (i+1):n_modes) {
            if(j <= n_modes && !modes_df$processed[j]) {
              if(abs(modes_df$Est_Mode[j] - current_mode) <= within) {
                # Add to same group
                modes_df$group[j] = group_counter
                modes_df$processed[j] = TRUE
              }
            }
          }
        }
        
        group_counter = group_counter + 1
      }
    }
    
    # Calculate group means
    result = modes_df %>%
      dplyr::filter(!is.na(group)) %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(Est_Mode = mean(Est_Mode), .groups = "drop") %>%
      dplyr::arrange(Est_Mode) %>%
      dplyr::mutate(group = 1:dplyr::n())  
    
  } else {
    # Original simple grouping
    df = df %>%
      dplyr::arrange(Est_Mode) %>%
      dplyr::mutate(group = cumsum(c(1, diff(Est_Mode) > within)))
    
    result = df %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(Est_Mode = mean(Est_Mode), .groups = "drop") %>%
      dplyr::arrange(Est_Mode)
  }
  
  return(result)
}

#' Fast MCMC for mixture models (Metropolis-Hastings-within-partial-Gibbs)
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
  result$method = "metropolis-hastings"
  
  return(result)
}

#' Dirichlet MCMC (identical interface to MM_MH)
#' 
#' @param y Numeric vector of data
#' @param grp Number of mixture components
#' @param prior_means Prior means for components
#' @param ids Vector of IDs for validation
#' @param n_iter Number of MCMC iterations
#' @param burnin Burn-in period
#' @param proposal_sd Proposal standard deviation
#' @param dirichlet_alpha Dirichlet concentration parameter
#' @param seed Random seed
#' @return List with MCMC results (SAME FORMAT as MM_MH)
#' @export
MM_MH_dirichlet <- function(y, grp, prior_means = NULL, ids,
                            n_iter = 5000, burnin = 2000,
                            proposal_sd = 0.15, dirichlet_alpha = 2.0,
                            seed = NULL) {
  
  # Same validation as MM_MH
  if(missing(ids)) {
    stop("'ids' argument is required for validation.")
  }
  
  if(is.null(ids)) {
    stop("'ids' cannot be NULL. Provide a vector of IDs for validation.")
  }
  
  if(length(ids) != length(y)) {
    stop("ID length (", length(ids), ") doesn't match y length (", length(y), ").")
  }
  
  if(is.null(seed)) {
    seed = sample.int(.Machine$integer.max, 1)
  }
  
  if(is.null(prior_means)) {
    prior_means = seq(min(y), max(y), length.out = grp)
  }
  
  if(length(prior_means) != grp) {
    warning("prior_means length doesn't match grp. Using equally spaced means.")
    prior_means = seq(min(y), max(y), length.out = grp)
  }
  
  # Call Dirichlet C++ function
  result = MM_MH_dirichlet_cpp(
    y = as.numeric(y),
    prior_means = as.numeric(prior_means),
    n_iter = as.integer(n_iter),
    burnin = as.integer(burnin),
    proposal_sd = as.numeric(proposal_sd),
    dirichlet_alpha = as.numeric(dirichlet_alpha),
    seed = as.integer(seed)
  )
  
  # Ensure compatibility with create_MM_output
  if(is.null(result$mean_acceptance)) {
    result$mean_acceptance = rep(0.5, result$n_components)
  }
  
  # Store IDs 
  result$ids = ids
  result$method = "dirichlet"  
  
  return(result)
}



#' Create output data frame
#' 
#' Converts MCMC results to exact CSV format
#' 
#' @param mcmc_result Output from MM_MH() or MM_MH_dirichlet()
#' @param y_original Original y values (if different from mcmc_result$y)
#' @param group_original Original group assignments (optional)
#' @param main_class Category/class name
#' @param max_groups Maximum number of groups for output columns
#' @return Data frame in CSV format
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

#' Parallel Bayesian mixture modeling using Markov Chain Monte Carlo (MCMC)
#' 
#' Performs multimodal probability assignment using either:
#' 1. Metropolis-Hastings-within-partial-Gibbs 
#' 2. Dirichlet-Multinomial 
#'  
#' @param data Input data frame
#' @param varCLASS Character, category variable name (required)
#' @param varY Character, value variable name (required)
#' @param varID Character, ID variable name (required)
#' @param method Density estimator method ("sj-dpi", "bcv", "ucv", "nrd") (default: "sj-dpi")
#' @param within Range parameter for grouping modes (default: 1.0)
#' @param maxNGROUP Maximum number of groups (default: 5)
#' @param out_dir Output directory for CSV files (if NULL, returns combined data frame)
#' @param n_workers Number of parallel workers (default: 3)
#' @param n_iter Number of MCMC iterations (default: 1000 for metropolis, 3000 for dirichlet)
#' @param burnin Burn-in period (default: 500 for metropolis, 1000 for dirichlet)
#' @param proposal_sd Proposal standard deviation for component means (default: 0.15)
#' @param sj_adjust Adjustment factor for bandwidth methods (default: 0.5, smaller -> more modes, higher -> fewer modes)
#' @param mcmc_method "metropolis" or "dirichlet"(default: "metropolis")
#' @param dirichlet_alpha Dirichlet concentration parameter (default: 2.0)
#' @return Data frame (if out_dir is NULL) or writes CSV files
#' @export
fuss_PARALLEL_mcmc <- function(data, 
                               varCLASS, 
                               varY, 
                               varID, 
                               method = "sj-dpi", 
                               within = 1.0, 
                               maxNGROUP = 5, 
                               out_dir = NULL, 
                               n_workers = 3,  
                               n_iter = NULL, 
                               burnin = NULL,
                               proposal_sd = 0.15, 
                               sj_adjust = 0.5,
                               mcmc_method = "metropolis",
                               dirichlet_alpha = 2.0) {
  
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
  
  # Set sensible defaults based on MCMC method
  if(is.null(n_iter)) {
    n_iter = ifelse(mcmc_method == "dirichlet", 3000, 1000)
  }
  
  if(is.null(burnin)) {
    burnin = ifelse(mcmc_method == "dirichlet", 1000, 500)
  }
  
  # Split data by category
  categories = unique(data[[varCLASS]])
  data_list = split(data, data[[varCLASS]])
  
  message("Processing ", length(categories), " categories in parallel with ", 
          n_workers, " workers")
  message("Settings: MCMC = ", mcmc_method, 
          ", Mode detection = ", method, 
          ", sj_adjust = ", sj_adjust,
          ", within = ", within)
  
  # Define processing function
  process_category <- function(cat_data, varCLASS, varY, varID,
                               method, within, maxNGROUP, out_dir,
                               n_iter, burnin, proposal_sd, sj_adjust,
                               mcmc_method, dirichlet_alpha) {
    
    cat_name = as.character(unique(cat_data[[varCLASS]])[1])
    message("Processing category: ", cat_name, " (", mcmc_method, ")")
    
    y = cat_data[[varY]]
    ids = cat_data[[varID]]
    
    # Skip if insufficient data
    if(length(y) < 5) {
      warning("  Skipping category '", cat_name, "': insufficient data (n = ", length(y), ")")
      return(NULL)
    }
    
    # Mode detection with enhanced method
    mode_result = tryCatch({
      get_MODES_enhanced(y, adjust = sj_adjust, threshold = 1.0)
    }, error = function(e) {
      warning("  Mode detection failed for '", cat_name, "': ", e$message)
      return(NULL)
    })
    
    if(is.null(mode_result)) {
      return(NULL)
    }
    
    if(method %in% names(mode_result)) {
      modes_df = mode_result[[method]]
    } else {
      modes_df = mode_result[["sj-dpi"]]
      warning("Method '", method, "' not found. Using 'sj-dpi'.")
    }
    
    # Group modes
    modes_grouped = group_MODES_enhanced(modes_df, within = within)
    n_components = nrow(modes_grouped)
    
    # Fallback if no modes detected
    if(n_components == 0) {
      warning("  No modes detected for '", cat_name, "'. Using 2 equally spaced modes.")
      n_components = 2
      prior_means = seq(min(y), max(y), length.out = n_components)
    } else {
      prior_means = modes_grouped$Est_Mode
    }
    
    message("  Detected ", n_components, " components")
    
    # Choose MCMC method
    mh_result = tryCatch({
      if(mcmc_method == "dirichlet") {
        MM_MH_dirichlet(
          y = y,
          grp = n_components,
          prior_means = prior_means,
          ids = ids,
          n_iter = n_iter,
          burnin = burnin,
          proposal_sd = proposal_sd,
          dirichlet_alpha = dirichlet_alpha,
          seed = 123
        )
      } else {
        MM_MH(
          y = y,
          grp = n_components,
          prior_means = prior_means,
          ids = ids,
          n_iter = n_iter,
          burnin = burnin,
          proposal_sd = proposal_sd,
          seed = 123
        )
      }
    }, error = function(e) {
      warning("  MCMC failed for '", cat_name, "': ", e$message)
      return(NULL)
    })
    
    if(is.null(mh_result)) {
      return(NULL)
    }
    
    message("  MCMC complete for '", cat_name, "'")
    
    # Create output
    output_df = create_MM_output(
      mcmc_result = mh_result,
      y_original = y,
      group_original = NULL,
      main_class = as.character(cat_name),
      max_groups = maxNGROUP
    )
    
    # Add method info
    attr(output_df, "mcmc_method") = mcmc_method
    
    # Write to CSV if out_dir provided
    if(!is.null(out_dir)) {
      if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
      filename = paste0("df_", gsub("[^a-zA-Z0-9]", "_", cat_name), "_", mcmc_method, ".csv")
      filepath = file.path(out_dir, filename)
      write.csv(output_df, filepath, row.names = FALSE, quote = TRUE)
      message("  Written: ", filepath)
      return(NULL)
    }
    
    return(output_df)
  }
  
  # Setup parallel processing
  future::plan(future::multisession, workers = n_workers)
  
  # Process each category in parallel
  result_list = furrr::future_map(data_list, function(cat_data) {
    
    process_category(
      cat_data = cat_data,
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
      sj_adjust = sj_adjust,
      mcmc_method = mcmc_method,
      dirichlet_alpha = dirichlet_alpha
    )
    
  }, .options = furrr::furrr_options(
    packages = "MultiModalR",
    seed = TRUE
  ), .progress = TRUE)
  
  # Reset to sequential
  future::plan(future::sequential)
  
  # If out_dir is NULL, combine and return results
  if(is.null(out_dir)) {
    result_list = result_list[!sapply(result_list, is.null)]
    
    if(length(result_list) > 0) {
      combined_result = do.call(rbind, result_list)
      message("Parallel analysis complete. Combined result has ", 
              nrow(combined_result), " rows.")
      
      # Add method attribute
      attr(combined_result, "mcmc_method") = mcmc_method
      attr(combined_result, "mode_method") = method
      attr(combined_result, "sj_adjust") = sj_adjust
      
      return(combined_result)
    } else {
      warning("No results were generated. Check your data and parameters.")
      return(NULL)
    }
  }
  
  message("Analysis complete. Results written to: ", out_dir)
  return(invisible(NULL))
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
  
  # ========== FIX: Calculate ALIGNED accuracy (handles label switching) ==========
  
  # For each category, align labels by value order and track correct counts
  all_correct = c()
  accuracy_by_class = list()
  
  for(category in unique(matched$Main_Class)) {
    cat_data = matched %>% dplyr::filter(Main_Class == category)
    
    # Get unique groups
    true_groups = sort(unique(cat_data[[subpop_col]]))
    pred_groups = sort(unique(cat_data$Assigned_Group))
    
    if(length(true_groups) != length(pred_groups)) {
      # If group counts don't match, use raw accuracy
      correct_vector = cat_data[[subpop_col]] == cat_data$Assigned_Group
    } else {
      # Calculate mean value for each true group
      true_means = cat_data %>%
        dplyr::group_by(!!rlang::sym(subpop_col)) %>%
        dplyr::summarize(mean_val = mean(!!rlang::sym(value_col)), .groups = 'drop') %>%
        dplyr::arrange(mean_val) %>%
        dplyr::pull(!!rlang::sym(subpop_col))
      
      # Calculate mean value for each predicted group  
      pred_means = cat_data %>%
        dplyr::group_by(Assigned_Group) %>%
        dplyr::summarize(mean_val = mean(y), .groups = 'drop') %>%
        dplyr::arrange(mean_val) %>%
        dplyr::pull(Assigned_Group)
      
      # Create mapping: align by value order
      mapping = setNames(true_means, pred_means)
      
      # Apply mapping to calculate aligned accuracy
      cat_data$Assigned_Aligned = mapping[as.character(cat_data$Assigned_Group)]
      correct_vector = cat_data[[subpop_col]] == cat_data$Assigned_Aligned
    }
    
    # Store correct/incorrect for overall calculation
    all_correct = c(all_correct, correct_vector)
    
    # Calculate accuracy for this category
    acc = round(mean(correct_vector, na.rm = TRUE) * 100, 1)
    
    accuracy_by_class[[category]] = data.frame(
      Main_Class = category,
      accuracy = acc,
      n = nrow(cat_data)
    )
  }
  
  # Calculate overall accuracy correctly
  overall_accuracy = round(mean(all_correct, na.rm = TRUE) * 100, 1)
  
  # Create label data for plot (using ALIGNED accuracy)
  label_data = purrr::map_dfr(accuracy_by_class, ~ .x) %>%
    dplyr::mutate(Percent = format(accuracy, nsmall = 1))
  
  # Create color palette for subgroups
  n_groups = length(unique(predicted$Assigned_Group))
  group_colors = c("firebrick2", "forestgreen", "cyan3", "gold", "purple", "orange")
  group_colors = group_colors[1:min(n_groups, length(group_colors))]
  
  # Plot validation - EXACT SAME STYLING AS BEFORE
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

