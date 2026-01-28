# R/create_dummy_data.R
#' Create multimodal dummy dataset
#' 
#' Generates the same dummy dataset used in the package. This is useful if users
#' want to generate similar data with different parameters.
#'
#' @param seed Random seed for reproducibility (default: 123)
#' @param n_categories Number of categories (default: 9)
#' @param n_per_group Number of observations per subgroup per category (default: 25)
#' @param n_subgroups Number of subgroups per category (default: 3)
#' @return A data frame with multimodal data
#' @export
#' @examples
#' # Generate the default dataset
#' df <- create_multimodal_dummy()
#' 
#' # Generate with different parameters
#' df2 <- create_multimodal_dummy(seed = 12345, n_categories = 6)
create_multimodal_dummy <- function(seed = 123, 
                                    n_categories = 9, 
                                    n_per_group = 25,
                                    n_subgroups = 3) {
  
  # Check if truncnorm is available
  if (!requireNamespace("truncnorm", quietly = TRUE)) {
    stop("Package 'truncnorm' is required to generate dummy data.
Please install it with: install.packages('truncnorm')")
  }
  
  set.seed(seed)
  
  # Create category names (AA, BB, CC, ...)
  categories <- rep(paste0(LETTERS[1:n_categories], LETTERS[1:n_categories]), 
                    each = n_per_group * n_subgroups)
  
  # Create subpopulation labels
  subpopulations <- rep(rep(paste("Group", 1:n_subgroups), each = n_per_group), 
                        times = n_categories)
  
  # Total observations
  n_total <- n_categories * n_subgroups * n_per_group
  
  # Generate values with varying means for each subgroup
  values <- numeric(n_total)
  
  for (i in 1:n_categories) {
    for (j in 1:n_subgroups) {
      # Calculate indices for this subgroup
      start_idx <- ((i - 1) * n_subgroups * n_per_group) + ((j - 1) * n_per_group) + 1
      end_idx <- start_idx + n_per_group - 1
      
      # Create a progression of means across subgroups
      base_mean <- 6 + (j - 1) * 1.5
      
      # Add slight variations across categories
      category_offset <- (i - 1) * 0.1
      subgroup_mean <- base_mean + category_offset
      
      # Generate truncated normal values
      values[start_idx:end_idx] <- truncnorm::rtruncnorm(
        n_per_group,
        a = 5,
        b = 10,
        mean = subgroup_mean,
        sd = 0.5 + (j - 1) * 0.1  # Increase variance for later groups
      )
    }
  }
  
  # Create data frame
  df <- data.frame(
    Category = factor(categories),
    Subpopulation = factor(subpopulations, levels = paste("Group", 1:n_subgroups)),
    Value = values,
    ID = 1:n_total
  )
  
  return(df)
}
