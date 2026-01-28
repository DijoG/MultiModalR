# data/multimodal_dummy.R
#' Multimodal Dummy Dataset for Testing MultiModalR
#'
#' A simulated dataset containing 9 categories (AA, BB, ..., II) each with 
#' 3 distinct subpopulations (Group 1, Group 2, Group 3) following 
#' truncated normal distributions.
#'
#' @format A data frame with 675 rows and 4 variables:
#' \describe{
#'   \item{Category}{Factor with 9 levels: AA, BB, CC, DD, EE, FF, GG, HH, II}
#'   \item{Subpopulation}{Factor with 3 levels: Group 1, Group 2, Group 3}
#'   \item{Value}{Numeric values between 5 and 10}
#'   \item{ID}{Unique identifier from 1 to 675}
#' }
#'
#' @details This dataset is useful for demonstrating the capabilities of 
#' MultiModalR package. Each category contains three distinct subpopulations
#' with overlapping but separable distributions, making it ideal for testing
#' multimodal mixture modeling algorithms.
#'
#' @source Generated with set.seed(5) using rtruncnorm package
#'
#' @examples
#' # Load the dataset
#' data(multimodal_dummy)
#' 
#' # View structure
#' str(multimodal_dummy)
#' 
#' # Summary statistics
#' summary(multimodal_dummy)
#' 
#' # Plot data
#' library(ggplot2)
#' ggplot(multimodal_dummy, aes(x = Value, fill = Subpopulation)) +
#'   geom_density(alpha = 0.5, color = NA) +
#'   facet_wrap(~Category) +
#'   theme_dark()
#' 
#' # Use with MultiModalR
#' \dontrun{
#' library(MultiModalR)
#' results <- fuss_PARALLEL_mcmc(
#'   data = multimodal_dummy,
#'   varCLASS = "Category",
#'   varY = "Value",
#'   varID = "ID",
#'   n_workers = 3
#' )
#' }
"multimodal_dummy"