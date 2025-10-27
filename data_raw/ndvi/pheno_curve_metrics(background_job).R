library(tidyverse)
library(phenofit)
library(oce)
library(progress)

# Load raw NDVI extraction data
read_rds("./raw_extracts/extracted_ndvi_value.rds") -> extracts

# # Extracting unique tree values
# extracts %>%
#   group_by(UID) %>%
#   summarise(n = n()) %>%
#   # sample_n(100) %>% # Randomly sample for testing
#   pull(UID) -> trees

read_rds("./extracted_metrics.rds") %>% 
  group_by(UID) %>% 
  summarise(n = n()) %>% 
  pull(UID) -> trees

# Initialize empty data frames with expected structure
pheno_metrix_date <- tibble()
pheno_metrix_doy <- tibble()
fit_stats <- tibble()

# Track successful and failed processing
success_count <- 0
fail_count <- 0

pb <- progress_bar$new(
  format = "Processing [:bar] :percent | Tree :current/:total | Elapsed: :elapsed | ETA: :eta",
  total = length(trees),
  clear = FALSE,
  width = 80
)

# Loop through each tree
for (i in trees) {
  
  tryCatch({
    
    # Filter data for current tree and preprocess NDVI
    extracts %>% 
      filter(UID == i) %>% 
      # Remove spikes/outliers using median-based despiking
      mutate(ndvi = despike(ndvi, "median", n = 3, k = 5)) %>%
      # Detrend NDVI
      mutate(
        time_numeric = as.numeric(im_date),
        trend = predict(lm(ndvi ~ time_numeric)),
        ndvi_detrended = ndvi - trend + mean(ndvi, na.rm = TRUE)
      ) -> extracts_samp
    
    # Prepare original NDVI data
    check_input(
      t = extracts_samp$im_date,
      y = extracts_samp$ndvi, 
      south = FALSE
    ) -> test_input
    
    # Prepare detrended NDVI data
    check_input(
      t = extracts_samp$im_date, 
      y = extracts_samp$ndvi_detrended,
      south = FALSE
    ) -> x_detrended 
    
    # Detect growing seasons - SET allow.cartesian=TRUE to fix join error
    season_mov(
      x_detrended,
      options = list(
        rFUN = "smooth_wHANTS",
        allow.cartesian = TRUE  # Allow cartesian joins
      )
    ) -> brks_mov 
    
    # Fit phenology curves
    fit <- curvefits(test_input, brks_mov)
    
    # Extract goodness of fit statistics
    stats <- get_GOF(fit)
    
    # Extract phenology metrics
    metrics <- get_pheno(fit, "Elmore")
    
    # Check if metrics were successfully extracted
    if (!is.null(metrics$date$Elmore) && !is.null(metrics$doy$Elmore)) {
      
      # Extract date-based metrics
      metrics_date <- metrics$date$Elmore %>% 
        mutate(UID = i)
      
      # Extract DOY-based metrics
      metrics_doy <- metrics$doy$Elmore %>% 
        mutate(UID = i)
      
      # Create lookup table with fit objects
      fit_lookup <- data.frame(
        flag = names(fit),
        fit_data = I(fit))
      
      # Add tree ID and join fit objects to statistics
      stats <- stats %>% 
        mutate(UID = i) %>% 
        left_join(fit_lookup, by = "flag")
      
      # Accumulate results
      pheno_metrix_date <- bind_rows(pheno_metrix_date, metrics_date)
      pheno_metrix_doy <- bind_rows(pheno_metrix_doy, metrics_doy)
      fit_stats <- bind_rows(fit_stats, stats)
      
      success_count <- success_count + 1
    }
    
    pb$tick()
    
  }, error = function(e) {
    message(paste("Skipping tree", i, "- Error:", e$message))
    fail_count <- fail_count + 1
    pb$tick()
  })
}

message(sprintf("\nProcessing complete: %d successful, %d failed", success_count, fail_count))

# Check if we have any successful results before proceeding
if (nrow(pheno_metrix_doy) == 0 || nrow(pheno_metrix_date) == 0) {
  stop("No trees were successfully processed. Cannot proceed with analysis.")
}

# Join DOY and date phenology metrics
phen_metrix <- left_join(
  pheno_metrix_doy, 
  pheno_metrix_date, 
  by = c("UID", "flag", "origin"),
  suffix = c("_doy", "_date"))

# Join fit statistics with phenology metrics
phen_metrix <- fit_stats %>% 
  filter(meth == "Elmore") %>%
  right_join(phen_metrix, by = c("UID", "flag"))

# Filter for high-quality trees
sample <- phen_metrix %>% 
  group_by(UID) %>% 
  summarise(
    mean_r2 = mean(R2, na.rm = TRUE),
    n = n()
  ) %>% 
  filter(n >= 6 & mean_r2 >= .9) %>%
  pull(UID)

message(sprintf("High-quality trees retained: %d out of %d successful", 
                length(sample), length(unique(phen_metrix$UID))))

# Filter and save
phen_metrix %>% 
  filter(UID %in% sample) %>% 
  write_rds("./pheno_curve_metrics/extracted_metrics.rds")