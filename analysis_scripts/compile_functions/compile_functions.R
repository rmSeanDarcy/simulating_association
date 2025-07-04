###########################################################################################
### SAMSARA - Loading Set of analyses data and preparing for plotting                   ###
###########################################################################################

###########################################################################################
# Load and prep data
compile_data <- function(experiment_name, results_filename = "results.csv") {
  # Build the experiment path from the current working directory
  experiment_path <- file.path(paste0(getwd(), '/simulation_data/', experiment_name))
  if (!dir.exists(experiment_path)) {
    stop("Experiment folder does not exist: ", experiment_path)
  }
  # Get full paths to all matching results files within the experiment folder
  result_paths <- list.files(
    path = experiment_path,
    pattern = results_filename,
    recursive = TRUE,
    full.names = TRUE
  )
  # Filter: only keep paths with at least 3 levels (e.g., treatment/simulation/results.csv)
  result_paths <- result_paths[sapply(strsplit(result_paths, .Platform$file.sep), length) >= 3]
  # Read and tag each results file
  result_list <- lapply(result_paths, function(file_path) {
    path_parts <- strsplit(file_path, .Platform$file.sep)[[1]]
    n <- length(path_parts)
    # Extract treatment and simulation folder names
    simulation <- path_parts[n - 1]
    treatment  <- path_parts[n - 2]
    df <- tryCatch(read.csv(file_path)[,-1], error = function(e) return(NULL))
    if (!is.null(df)) {
      df$treatment <- treatment
      df$simulation <- simulation
      df$experiment <- experiment_name  # Add experiment name too, if useful
    }
    return(df)
  })
  # Filter out failed loads
  result_list <- Filter(Negate(is.null), result_list)
  # Combine everything into one big data frame
  combined_df <- do.call(rbind, result_list)
  return(combined_df)
}




