#' Evaluate Phenology Method Combinations
#'
#' Tests combinations of smoothing functions and phenology extraction methods
#' on a random sample of trees to determine the proportion that meet quality
#' criteria (minimum observations and R² threshold).
#'
#' @param extracts A data frame containing NDVI time series data. Must include
#'   columns: UID (tree identifier), im_date (image date), and ndvi (NDVI values).
#' @param n_samples Integer specifying the number of trees to randomly sample
#'   for evaluation.
#' @param pheno_curves Character vector of phenology curve fitting methods to test.
#'   Options include: "AG", "Beck", "Elmore", "Zhang", "Gu", "Klos", "Doublelogistic.zhang".
#' @param season_curves Character vector of smoothing functions to test (without
#'   the "smooth_" prefix). Options include: "wHANTS", "wSG", "wWHIT", "wBisquare",
#'   "wChen", "wKong".
#' @param min_observations Integer specifying the minimum number of observations
#'   (seasons/years) required per tree. Default is 6.
#' @param min_r2 Numeric specifying the minimum mean R² threshold for including
#'   a tree in the "usable" count. Default is 0.9.
#'
#' @details
#' ## Processing Steps
#' For each combination of smoothing and phenology methods, the function:
#' \enumerate{
#'   \item Randomly samples \code{n_samples} trees from the dataset
#'   \item Applies despiking, detrending, and smoothing to NDVI time series
#'   \item Fits phenology curves using the specified method
#'   \item Calculates goodness-of-fit statistics
#'   \item Determines the proportion of trees meeting quality criteria
#' }
#' 
#' ## Quality Criteria
#' A tree is considered "usable" if it meets BOTH conditions:
#' \itemize{
#'   \item Has at least \code{min_observations} successful seasonal fits
#'   \item Has a mean R² across seasons >= \code{min_r2}
#' }
#' 
#' ## Smoothing Functions
#' The function tests different smoothing approaches:
#' \itemize{
#'   \item \strong{wHANTS}: Harmonic ANalysis of Time Series (weighted)
#'   \item \strong{wSG}: Savitzky-Golay filter (weighted)
#'   \item \strong{wWHIT}: Whittaker smoother (weighted)
#' }
#' 
#' ## Phenology Methods
#' Different curve fitting approaches:
#' \itemize{
#'   \item \strong{Elmore}: Double logistic with background trend
#'   \item \strong{AG}: Asymmetric Gaussian
#'   \item \strong{Beck}: Double logistic
#'   \item \strong{Zhang}: Asymmetric logistic
#' }
#'
#' @return A numeric matrix with smoothing methods as rows and phenology methods
#'   as columns. Each cell contains the proportion (0-1) of sampled trees that
#'   met the quality criteria for that method combination.
#'
#' @examples
#' # Define methods to test
#' pheno_curves <- c("AG", "Beck", "Elmore", "Zhang")
#' season_curves <- c("wHANTS", "wSG", "wWHIT")
#' 
#' # Test with 100 random trees
#' usability_matrix <- evaluate_pheno_methods(
#'   extracts = extracts,
#'   n_samples = 100,
#'   pheno_curves = pheno_curves,
#'   season_curves = season_curves,
#'   min_observations = 6,
#'   min_r2 = 0.9
#' )
#' 
#' # View results
#' print(usability_matrix)
#' 
#' # Test with stricter criteria
#' usability_strict <- evaluate_pheno_methods(
#'   extracts = extracts,
#'   n_samples = 500,
#'   pheno_curves = pheno_curves,
#'   season_curves = season_curves,
#'   min_observations = 8,
#'   min_r2 = 0.95
#' )
#' 
#' # Format results as percentages
#' usability_pct <- apply(usability_matrix, c(1, 2), 
#'                       function(x) sprintf("%.1f%%", x * 100))
#' print(usability_pct)
#'
#' @note
#' This function can be time-consuming for large sample sizes as it fits
#' multiple models for each tree. Progress messages are printed to track
#' completion. Trees that produce errors during processing are silently
#' skipped and not counted in the results.
#'
#' @seealso 
#' \code{\link{filter_pheno_data}} for processing all trees with a specific method
#' \code{\link{plot_year}} for visualizing individual tree phenology
#'
#' @references
#' Kong, D., et al. (2020). A robust method for reconstructing global MODIS EVI
#' time series on the Google Earth Engine. \emph{ISPRS Journal of Photogrammetry
#' and Remote Sensing}, 155, 13-24.
#'
#' @export
evaluate_pheno_methods <- function(extracts, n_samples, pheno_curves, season_curves, 
                                   min_observations = 6, min_r2 = 0.9) {
  
  # Random sample of trees
  trees <- extracts %>% 
    group_by(UID) %>% 
    summarise(n = n(), .groups = "drop") %>% 
    sample_n(n_samples) %>% 
    pull(UID)
  
  message(paste("Randomly sampled", n_samples, "trees for analysis\n"))
  
  # Initialize results matrix
  results <- matrix(NA, 
                    nrow = length(season_curves), 
                    ncol = length(pheno_curves),
                    dimnames = list(season_curves, pheno_curves))
  
  total_trees <- length(trees)
  
  # Loop through each combination
  for (sc in season_curves) {
    for (pc in pheno_curves) {
      
      message(paste("\nProcessing:", sc, "/", pc))
      
      # Initialize data frames for this combination
      pheno_metrix_date <- data.frame()
      pheno_metrix_doy <- data.frame()
      fit_stats <- data.frame()
      
      # Loop through trees
      for (i in trees) {
        tryCatch({
          
          extracts %>% 
            filter(UID == i) %>% 
            mutate(ndvi = despike(ndvi, "median", n = 3, k = 5)) %>%
            mutate(
              time_numeric = as.numeric(im_date),
              trend = predict(lm(ndvi ~ time_numeric)),
              ndvi_detrended = ndvi - trend + mean(ndvi, na.rm = TRUE)
            ) -> extracts_samp
          
          check_input(
            t = extracts_samp$im_date,
            y = extracts_samp$ndvi, 
            south = FALSE
          ) -> test_input
          
          check_input(
            t = extracts_samp$im_date, 
            y = extracts_samp$ndvi_detrended,
            south = FALSE
          ) -> x_detrended 
          
          season_mov(
            x_detrended,
            options = list(
              rFUN = paste0("smooth_", sc)
            )
          ) -> brks_mov 
          
          curvefits(test_input, brks_mov) -> fit
          
          get_GOF(fit) -> stats
          get_pheno(fit, pc) -> metrics
          
          # Extract date and doy metrics
          metrics$date[[pc]] %>% 
            mutate(UID = i) -> metrics_date
          
          metrics$doy[[pc]] %>% 
            mutate(UID = i) -> metrics_doy
          
          stats %>% 
            mutate(UID = i) -> stats
          
          # Bind rows
          pheno_metrix_date <- bind_rows(pheno_metrix_date, metrics_date)
          pheno_metrix_doy <- bind_rows(pheno_metrix_doy, metrics_doy)
          fit_stats <- bind_rows(fit_stats, stats)
          
        }, error = function(e) {
          # Silent error handling
        })
      }
      
      # Calculate proportion meeting criteria
      if (nrow(fit_stats) > 0) {
        
        phen_metrix <- left_join(
          pheno_metrix_doy, 
          pheno_metrix_date, 
          by = c("UID", "flag", "origin"),
          suffix = c("_doy", "_date")
        )
        
        fit_stats_filtered <- fit_stats %>% 
          filter(meth == pc) %>% 
          group_by(UID) %>% 
          summarise(mean_r2 = mean(R2, na.rm = TRUE), .groups = "drop")
        
        n_meeting_criteria <- phen_metrix %>% 
          count(UID) %>% 
          left_join(fit_stats_filtered, by = "UID") %>% 
          filter(n >= min_observations & mean_r2 >= min_r2) %>% 
          nrow()
        
        results[sc, pc] <- n_meeting_criteria / total_trees
        
        message(paste(sc, "/", pc, ":", 
                      sprintf("%.2f", results[sc, pc]), "usable (",
                      n_meeting_criteria, "of", total_trees, "trees)"))
      } else {
        results[sc, pc] <- 0
        message(paste(sc, "/", pc, ": 0.00 usable (no successful fits)"))
      }
    }
  }
  
  return(results)
}