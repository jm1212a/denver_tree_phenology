library(tidyverse)
library(pracma)

plot_growing_season <- function(ndvi_data, ndvi_obs , df = 30, maxtroughhight = .3, 
                                mintroughdistance = 300, max_season_length = 380){
  
  ndvi_data %>% 
    pull({{ ndvi_obs }}) -> y_vals
  
  smooth.spline(as.numeric(ndvi_data$im_date), 
                y_vals, 
                df = df) -> sm 
  
  seq(min(as.numeric(ndvi_data$im_date)), 
      max(as.numeric(ndvi_data$im_date)), 
      by = 1) -> x_daily 
  
  predict(sm, x_daily)$y -> y_smooth 
  
  findpeaks(
    -y_smooth,
    minpeakdistance = mintroughdistance,
    minpeakheight = -maxtroughhight) -> troughs
  
  sort(as.Date(x_daily[troughs[, 2]], origin = "1970-01-01")) -> trough_dates 
  
  c(min(ndvi_data$im_date), trough_dates) -> trough_dates 
  
  plot(extracts_samp$im_date, extracts_samp$ndvi)
  
  lines(sm, col = "red")
  lines(as.Date(x_daily, origin = "1970-01-01"), y_smooth, col = "red")
  abline(v = as.numeric(trough_dates), col = "blue", lty = 2)
  
}

# plot_growing_season(extracts_samp, ndvi_obs = ndvi, maxtroughhight = .4) #use example 