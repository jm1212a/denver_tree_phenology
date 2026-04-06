library(tidyverse)
library(pracma)
library(minpack.lm) 

source("../../R/growing_season.R")
source("../../R/fit_annual_curve.R")

pheno_metrics <- function(location, ndvi_data) {
  ndvi_data %>%
    filter(UID == location) %>% 
    mutate(ndvi_normalized = despike(ndvi_normalized, "median", n = 2, k = 5)) -> extracts_samp
  
  growing_season(extracts_samp, ndvi_obs = ndvi_normalized, maxtroughhight = .6) -> seasons
  
  map(seasons, possibly(~fit_annual_curve(.x, ndvi_obs = ndvi_normalized), otherwise = NULL)) %>%
    compact() %>%
    bind_rows()
}