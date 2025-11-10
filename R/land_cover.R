library(tidyverse)
library(sf)
library(tigris)
library(gridExtra)
library(rstac)
library(terra)
library(tidyterra)

stac_source <- stac("https://planetarycomputer.microsoft.com/api/stac/v1/")

denver_county <- counties(state = "CO", cb = TRUE) %>%
  filter(NAME == "Denver") %>%
  st_transform(crs = 4326) # adjust to local projection for final sample

stac_query <- stac_search(
  q = stac_source,
  collections = "drcog-lulc",
  intersects = denver_county
)

signed_stac_query <- items_sign(
  post_request(stac_query),
  sign_planetary_computer()
)

class_df <- map_dfr(signed_stac_query$features[[1]]$assets$data$`classification:classes`, as_tibble)

class_df$col <- paste0("#", class_df$color_hint)

ref_DLCD <- vector("list", length = 93)  # preallocate list

for (i in 1:93) {
  # Load raster
  r <- rast(paste0("/vsicurl/", signed_stac_query$features[[i]]$assets$data$href))
  
  # Get raster corners
  e <- ext(r)
  corners <- vect(rbind(
    c(xmin(e), ymin(e)),
    c(xmin(e), ymax(e)),
    c(xmax(e), ymin(e)),
    c(xmax(e), ymax(e))
  ), crs = crs(r))
  
  corners_ll <- project(corners, "EPSG:4326")
  
  # Store metadata for this raster
  ref_DLCD[[i]] <- data.frame(
    rast_id = paste0("/vsicurl/", signed_stac_query$features[[i]]$assets$data$href),
    xmin = xmin(corners_ll),
    xmax = xmax(corners_ll),
    ymin = ymin(corners_ll),
    ymax = ymax(corners_ll)
  )
}

# Bind list into one data frame
ref_DLCD <- dplyr::bind_rows(ref_DLCD)

land_cover <- function(geom) {
  
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