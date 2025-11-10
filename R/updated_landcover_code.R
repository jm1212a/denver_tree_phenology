library(tidyverse)
library(sf)
library(tigris)
library(gridExtra)
library(rstac)
library(terra)
library(tidyterra)

# Land cover function with mosaicking for edge trees
# Returns ONLY numeric values - NO LISTS
land_cover_mosaic <- function(geom) {
  
  # wrap raw geometry in sfc with CRS
  if (!inherits(geom, "sfc")) {
    geom <- st_sfc(geom, crs = 4326)
  }
  
  coords <- st_coordinates(geom)[1, ]
  x <- coords[1]; y <- coords[2]
  
  # Buffer extent in degrees (90m buffer)
  buffer_deg <- 90 / 111320
  
  # Calculate bounding box of buffered point
  buffer_bbox <- data.frame(
    xmin = x - buffer_deg,
    xmax = x + buffer_deg,
    ymin = y - buffer_deg,
    ymax = y + buffer_deg
  )
  
  # Find all rasters that intersect with the buffered point
  intersecting_rasters <- ref_DLCD %>%
    filter(
      xmax >= buffer_bbox$xmin &
        xmin <= buffer_bbox$xmax &
        ymax >= buffer_bbox$ymin &
        ymin <= buffer_bbox$ymax
    )
  
  # Load raster(s)
  if (nrow(intersecting_rasters) == 0) {
    stop("No raster found for this point")
  } else if (nrow(intersecting_rasters) == 1) {
    r <- rast(intersecting_rasters$rast_id[1])
  } else {
    # Multiple rasters - load and mosaic them
    raster_list <- lapply(intersecting_rasters$rast_id, rast)
    
    ref_raster <- raster_list[[1]]
    
    if (length(raster_list) > 1) {
      for (i in 2:length(raster_list)) {
        raster_list[[i]] <- resample(raster_list[[i]], ref_raster, method = "near")
      }
    }
    
    r <- do.call(mosaic, raster_list)
  }
  
  # Extract and process
  bufs <- buffer(vect(geom), width = 90, quadsegs = 4)
  cell_counts <- terra::extract(r, bufs, touches = FALSE)
  
  # Ensure it's a proper data frame
  cell_counts <- as.data.frame(cell_counts)
  
  # Get column name for land cover values
  lc_col <- names(cell_counts)[2]
  
  # Group and count
  cell_counts <- cell_counts %>% 
    group_by(ID, !!sym(lc_col)) %>%
    summarise(count = n(), .groups = "drop")
  
  names(cell_counts)[2] <- "value"
  
  # Join with class definitions
  cell_counts <- cell_counts %>% 
    left_join(class_df, by = "value") %>% 
    select(ID, value, description, count) %>%
    mutate(description = if_else(is.na(description), "NA", as.character(description)))
  
  # Calculate total
  total_count <- sum(cell_counts$count)
  
  # Aggregate by description to prevent duplicates
  cell_agg <- cell_counts %>%
    group_by(description) %>%
    summarise(
      count = sum(count),
      per_surface = sum(count) / total_count,
      .groups = "drop"
    )
  
  # Build result manually - NO PIVOT_WIDER
  result <- tibble(ID = 1)
  
  for (i in 1:nrow(cell_agg)) {
    desc <- cell_agg$description[i]
    result[[paste0("%_of_", desc)]] <- as.numeric(cell_agg$per_surface[i])
    result[[paste0("#_of_pix_", desc)]] <- as.numeric(cell_agg$count[i])
  }
  
  return(result)
}
