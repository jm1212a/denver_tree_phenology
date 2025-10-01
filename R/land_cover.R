library(tidyverse)
library(sf)
library(tigris)
library(gridExtra)
library(rstac)
library(terra)
library(tidyterra)

stac_query <- stac_search(
  q = stac("https://planetarycomputer.microsoft.com/api/stac/v1/"),
  collections = "drcog-lulc"
)

land_cover <- function(geom) {
  
  ref_DLCD <- read_csv("../data_raw/mis._data_sets/ref_DLCD.csv")
  class_df <- read_csv("../data_raw/mis._data_sets/class_df.csv")
  
  # wrap raw geometry in sfc with CRS
  if (!inherits(geom, "sfc")) {
    geom <- st_sfc(geom, crs = 4326)
  }
  
  rast_id_match <- NA
  coords <- st_coordinates(geom)[1, ]
  x <- coords[1]; y <- coords[2]
  
  for (i in seq_len(nrow(ref_DLCD))) {
    bbox <- ref_DLCD[i, ]
    if (x >= bbox$xmin && x <= bbox$xmax &&
        y >= bbox$ymin && y <= bbox$ymax) {
      rast_id_match <- bbox$rast_id
      break
    }
  }
  
  r <- rast(rast_id_match)
  
  bufs <- buffer(vect(geom), width = 90, quadsegs = 4)
  cell_counts <- terra::extract(r, bufs)
  
  cell_counts %>% 
    bind_rows() -> cell_counts 
  
  cell_counts %>% 
    group_by(ID, cell_counts[,2]) %>%
    summarise(count = n(), .groups = "drop") -> cell_counts
  
  names(cell_counts)[2] <- "value"
  
  cell_counts %>% 
    left_join(class_df) %>% 
    select(ID, value,description,count) %>% 
    left_join(cell_counts %>% 
                group_by(ID) %>% 
                summarise(total = sum(count), .groups = "drop")) -> cell_counts
  
  cell_counts %>% 
    mutate(per_surface = count/total) %>% 
    select(-value, -count, -total) %>% 
    mutate(description_per = paste0("%_of_",description))%>% 
    select(-description) %>% 
    pivot_wider(names_from = description_per,values_from = per_surface, 
                values_fill = 0) -> cell_per
  
  cell_counts %>% 
    select(-value, -total) %>% 
    mutate(description_count = paste0("#_of_pix_",description)) %>% 
    select(-description) %>% 
    pivot_wider(names_from = description_count, values_from = count, 
                values_fill = 0) -> cell_count
  
  left_join(cell_per, cell_count) -> cell_summary
  
  return(cell_summary)
}