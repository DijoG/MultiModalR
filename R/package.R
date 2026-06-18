# Package-level documentation and imports
# This file is used by roxygen2 to generate the NAMESPACE file

#' @importFrom stats density setNames
#' @importFrom utils install.packages installed.packages write.csv
#' @importFrom dplyr %>%
NULL

# Declare global variables used in dplyr/data.table pipelines
# This prevents "no visible binding for global variable" notes
utils::globalVariables(c(
  "Density_Height",
  "group",
  "Est_Mode",
  "Main_Class",
  "mean_val",
  "Assigned_Group",
  "y",
  "accuracy",
  "Percent"
))
