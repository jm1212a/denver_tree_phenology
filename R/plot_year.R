library(tidyverse)

#' Plot Phenology Metrics for a Single Tree and Year
#'
#' Creates a scatter plot of NDVI observations with fitted Elmore curve and 
#' optional phenology metric markers for a specific tree and year combination.
#'
#' @param data A data frame containing phenology metrics and fit data. Must include
#'   columns: flag, UID, fit_data (list column), and phenology metric columns 
#'   (e.g., Greenup_date, Greenup_doy, etc.)
#' @param year Character string specifying the flag/year to plot (e.g., "2018_1")
#' @param tree_id Character or numeric specifying the tree UID to plot
#' @param metrics Character vector of metric names or metric groups to display.
#'   Default is an empty vector (no metrics displayed). See Details for available
#'   metric groups.
#'
#' @details
#' ## Metric Groups
#' The function recognizes the following metric groups:
#' \itemize{
#'   \item \strong{Inflection}: Greenup, Maturity, Senescence, Dormancy
#'   \item \strong{Derivative}: DER.sos, DER.pos, DER.eos
#'   \item \strong{GU}: UD, SD, DD, RD (curvature-based metrics)
#'   \item \strong{TRS2}: TRS2.sos, TRS2.eos (20% threshold)
#'   \item \strong{TRS5}: TRS5.sos, TRS5.eos (50% threshold)
#'   \item \strong{TRS6}: TRS6.sos, TRS6.eos (60% threshold)
#'   \item \strong{TRS}: All threshold-based metrics (TRS2, TRS5, TRS6)
#' }
#' 
#' Individual metric names can also be specified directly (e.g., "Greenup", "DER.sos").
#' 
#' ## Color Scheme
#' Each metric type has a predefined color for consistency across plots:
#' \itemize{
#'   \item Inflection metrics: green to orange to brown gradient
#'   \item Derivative metrics: pink shades
#'   \item GU metrics: yellow/gold shades
#'   \item TRS metrics: blue, cyan, and purple for different thresholds
#' }
#'
#' @return Invisibly returns NULL. The function is called for its side effect of
#'   creating a plot in the active graphics device.
#'
#' @examples
#' # Plot just the fitted curve (no metrics)
#' plot_year(phen_metrix, "2018_1", "20788")
#' 
#' # Plot with inflection point metrics
#' plot_year(phen_metrix, "2018_1", "20788", metrics = "Inflection")
#' 
#' # Plot with derivative-based metrics
#' plot_year(phen_metrix, "2018_1", "20788", metrics = "Derivative")
#' 
#' # Plot with multiple metric groups
#' plot_year(phen_metrix, "2018_1", "20788", 
#'           metrics = c("Inflection", "Derivative"))
#' 
#' # Plot all TRS thresholds
#' plot_year(phen_metrix, "2018_1", "20788", metrics = "TRS")
#' 
#' # Mix groups and individual metrics
#' plot_year(phen_metrix, "2018_1", "20788", 
#'           metrics = c("Inflection", "DER.sos", "TRS5"))
#' 
#' # Plot everything
#' plot_year(phen_metrix, "2018_1", "20788", 
#'           metrics = c("Inflection", "Derivative", "GU", "TRS"))
#'
#' @seealso 
#' \code{\link{filter_pheno_data}} for filtering phenology data by quality metrics
#'
#' @export
plot_year <- function(data, year, tree_id, 
                      metrics = c()){
  
  # Define metric groups
  metric_groups <- list(
    "Inflection" = c("Greenup", "Maturity", "Senescence", "Dormancy"),
    "Derivative" = c("DER.sos", "DER.pos", "DER.eos"),
    "GU" = c("UD", "SD", "DD", "RD"),
    "TRS2" = c("TRS2.sos", "TRS2.eos"),
    "TRS5" = c("TRS5.sos", "TRS5.eos"),
    "TRS6" = c("TRS6.sos", "TRS6.eos"),
    "TRS" = c("TRS2.sos", "TRS2.eos", "TRS5.sos", "TRS5.eos", "TRS6.sos", "TRS6.eos")
  )
  
  # Expand groups to individual metrics
  expanded_metrics <- c()
  for(m in metrics) {
    if(m %in% names(metric_groups)) {
      expanded_metrics <- c(expanded_metrics, metric_groups[[m]])
    } else {
      expanded_metrics <- c(expanded_metrics, m)
    }
  }
  
  data %>% 
    filter(.data$year == !!year) %>% 
    filter(UID == tree_id) -> data_year
  
  data_year %>% 
    pull(fit_data) %>% 
    .[[1]] -> plot_list
  
  elmore_function <- function(t, mn, mx, sos, rsp, eos, rau, m7) {
    mn + (mx - m7 * t) * (1/(1 + exp(-rsp * (t - sos))) - 
                            1/(1 + exp(-rau * (t - eos))))
  }
  
  as.Date(plot_list$data$t, origin = "2000-01-01") -> t_dates
  
  plot(t_dates, plot_list$data$y,
       main = paste("Tree:", tree_id, "Year:", year),
       xlab = "Date", ylab = "NDVI",
       pch = 16, col = "gray60",
       ylim = c(min(plot_list$data$y, na.rm = TRUE), 
                max(plot_list$data$y, na.rm = TRUE) * 1.15))
  
  lines(t_dates, 
        elmore_function(plot_list$data$t, 
                        plot_list$model$Elmore$par[1], 
                        plot_list$model$Elmore$par[2], 
                        plot_list$model$Elmore$par[3],
                        plot_list$model$Elmore$par[4], 
                        plot_list$model$Elmore$par[5], 
                        plot_list$model$Elmore$par[6], 
                        plot_list$model$Elmore$par[7]),
        col = "red", lwd = 2)
  
  # Define colors for all metrics
  colors <- c(
    "Greenup" = "lightgreen", 
    "Maturity" = "darkgreen", 
    "Senescence" = "orange", 
    "Dormancy" = "brown",
    "TRS2.sos" = "blue",
    "TRS2.eos" = "darkblue",
    "TRS5.sos" = "cyan",
    "TRS5.eos" = "darkcyan",
    "TRS6.sos" = "purple",
    "TRS6.eos" = "purple4",
    "DER.sos" = "pink",
    "DER.pos" = "hotpink",
    "DER.eos" = "deeppink",
    "UD" = "yellow3",
    "SD" = "gold",
    "DD" = "goldenrod",
    "RD" = "darkgoldenrod"
  )
  
  # Plot each metric
  y_offset <- 0.95
  for(metric in expanded_metrics) {
    date_col <- paste0(metric, "_date")
    doy_col <- paste0(metric, "_doy")
    
    if(date_col %in% names(data_year)) {
      metric_date <- data_year[[date_col]]
      metric_doy <- data_year[[doy_col]]
      metric_color <- colors[metric]
      if(is.na(metric_color)) metric_color <- "black"
      
      abline(v = metric_date, col = metric_color, lty = 2, lwd = 1.5)
      text(metric_date, par("usr")[4] * y_offset, 
           paste(metric, "\n", metric_date, "\nDOY:", metric_doy), 
           pos = 2, col = metric_color, cex = 0.7)
      
      y_offset <- y_offset - 0.05
    }
  }
  
  invisible(NULL)
}
