library(tidyverse)
library(pracma)


growing_season <- function(ndvi_data, ndvi_obs, df = 30, maxtroughhight = .3, 
                           mintroughdistance = 300, max_season_length = 380){
  
  ndvi_data %>% 
    pull({{ ndvi_obs }}) -> y_vals
  
  smooth.spline(as.numeric(ndvi_data$im_date), y_vals, df = df) -> sm 
  
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
  
  seasons <- list()
  
  for (i in seq_along(trough_dates[-length(trough_dates)])) {
    
    as.numeric(trough_dates[i+1] - trough_dates[i]) -> gap 
    
    if (gap < max_season_length) {
      trough_dates[i] -> season_start 
      trough_dates[i + 1] -> season_end
      
      ndvi_data %>%
        filter(im_date >= season_start & im_date < season_end) %>%
        mutate(season = year(season_start),
               flag = paste0(lubridate::year(season_start), "_", 1))%>% 
        mutate(base_date = date(str_c(season,"-","01","-","01")),
               im_doy = as.numeric(difftime(im_date, base_date, units = "days"))) -> seasons[[i]] 
    }
  }
  return(seasons)
}