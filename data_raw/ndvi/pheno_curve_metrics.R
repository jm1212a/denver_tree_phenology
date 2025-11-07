library(tidyverse)
library(phenofit)
library(oce)
library(zoo)
library(furrr)
library(progressr)

# Load raw NDVI extraction data
read_rds("./raw_extracts/extracted_ndvi_value.rds") -> extracts

# Sample trees for testing
extracts %>% 
  group_by(UID) %>% 
  summarise(n = n()) %>% 
  sample_n(100) %>%
  pull(UID) -> trees

# Pre-split data by tree 
extracts %>% 
  filter(UID %in% trees) %>% 
  group_split(UID, .keep = TRUE) -> extracts_split 

total_trees <- length(trees)
message(sprintf("Total trees to process: %d", total_trees))

# Function to process one tree
process_tree <- function(tree_data, p = NULL) {
  
  i <- unique(tree_data$UID)
  
  tryCatch({
    
    # Preprocess NDVI
    tree_data %>% 
      mutate(ndvi = despike(ndvi, "median", n = 3, k = 5)) %>%
      mutate(
        time_numeric = as.numeric(im_date),
        trend = predict(lm(ndvi ~ time_numeric)),
        ndvi_detrended = ndvi - trend + mean(ndvi, na.rm = TRUE)
      ) -> extracts_samp
    
    # Prepare inputs
    test_input <- check_input(
      t = extracts_samp$im_date,
      y = extracts_samp$ndvi, 
      south = FALSE
    )
    
    x_detrended <- check_input(
      t = extracts_samp$im_date, 
      y = extracts_samp$ndvi_detrended,
      south = FALSE
    )
    
    # Detect growing seasons
    brks_mov <- season_mov(
      x_detrended,
      options = list(rFUN = "smooth_wWHIT")
    )
    
    # Fit phenology curves
    fit <- curvefits(test_input, brks_mov)
    stats <- get_GOF(fit)
    
    # Check fit quality
    fit_check <- stats %>% 
      filter(meth == "Elmore") %>% 
      mutate(UID = i) %>%
      group_by(UID) %>% 
      summarise(
        mean_r2 = mean(R2, na.rm = TRUE),  
        n = n()
      )
    
    if (fit_check$mean_r2 > 0.90 & fit_check$n >= 6) {
      
      # Extract phenology metrics
      metrics <- get_pheno(fit, "Elmore")
      
      metrics_date <- metrics$date$Elmore %>% 
        mutate(UID = i)
      
      metrics_doy <- metrics$doy$Elmore %>% 
        mutate(UID = i)
      
      # Create lookup tables
      fit_lookup <- data.frame(
        flag = names(fit),
        fit_data = I(fit))
      
      season_curve <- tibble(
        UID = i,
        season_data = list(test_input))
      
      # Join fit objects to statistics
      stats <- stats %>% 
        mutate(UID = i) %>%
        left_join(season_curve, by = "UID") %>% 
        left_join(fit_lookup, by = "flag")
      
      if (!is.null(p)) p(message = sprintf("Tree %s: Success", i))
      
      # Return results
      return(list(
        metrics_date = metrics_date,
        metrics_doy = metrics_doy,
        stats = stats
      ))
      
    } else {
      if (!is.null(p)) p(message = sprintf("Tree %s: Skipped (r2 < .9 or n < 6)", i))
      return(list(status = "skipped", UID = i))
    }
    
  }, error = function(e) {
    if (!is.null(p)) p(message = sprintf("Tree %s: Error - %s", i, e$message))
    return(list(status = "error", UID = i, error = e$message))
  })
}

# Set up parallel processing (use all cores minus 1)
plan(multisession, workers = availableCores() - 1)

handlers(global = TRUE)
handlers("progress")

start_time <- Sys.time()

# Process all trees in parallel with progress bar
with_progress({
  p <- progressor(steps = length(extracts_split))
  
  results <- future_map(
    extracts_split, 
    ~process_tree(.x, p),
    .options = furrr_options(seed = TRUE)
  )
})

end_time <- Sys.time()
message(sprintf("\nProcessing complete! Time elapsed: %.2f minutes", 
                as.numeric(difftime(end_time, start_time, units = "mins"))))

# Analyze results
successful <- sum(map_chr(results, ~.x$status %||% "success") == "success")
skipped <- sum(map_chr(results, ~.x$status %||% "success") == "skipped")
errors <- sum(map_chr(results, ~.x$status %||% "success") == "error")

message(sprintf("Results: %d successful, %d skipped, %d errors", 
                successful, skipped, errors))

results_successful <- results[map_chr(results, ~.x$status %||% "success") == "success"]

pheno_metrix_date <- map_df(results_successful, "metrics_date")
pheno_metrix_doy <- map_df(results_successful, "metrics_doy")
fit_stats <- map_df(results_successful, "stats")

# Join DOY and date phenology metrics
left_join(
  pheno_metrix_doy, 
  pheno_metrix_date, 
  by = c("UID", "flag", "origin"),
  suffix = c("_doy", "_date")) -> phen_metrix 

# Join fit statistics with phenology metrics
fit_stats %>% 
  filter(meth == "Elmore") %>%
  right_join(phen_metrix) %>% 
  write_rds("./pheno_curve_metrics/extracted_metrics_test.rds")

