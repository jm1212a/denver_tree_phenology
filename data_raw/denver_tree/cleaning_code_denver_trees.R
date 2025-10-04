library(tidyverse)
library(sf)
library(tigris)

denver_county <- counties(state = "CO", cb = TRUE) %>%
  filter(NAME == "Denver") %>%
  st_transform(crs = 4326) # adjust to local projection for final sample

source("../../R/land_cover.R")

read_rds("tree_inventory_denver_20250911.rds") %>%
  # Fix formatting of common names: convert "last, first" to "first last"
  mutate(
    species_common = if_else(
      str_detect(species_common, ","),# check if a comma exists
      str_replace(species_common, "^(.*),\\s*(.*)$", "\\2 \\1"), # swap parts around
      species_common  # leave unchanged if no comma
    ),
    # Combine scientific and common names into one column
    botanic_common = paste(species_botanic, "(", species_common, ")"),
    # Standardize disease/pest field: anything containing "N/A" becomes exactly "N/A"
    disease_pest_def = if_else(str_detect(disease_pest_def, "N/A"), 
                               "N/A", disease_pest_def)
  ) %>%
  filter(species_botanic != "Tree requested") %>% # Remove rows with placeholder species
  # Extract the Genus + species (binomial name) from species_botanic
  mutate(
    species_sci = str_extract(species_botanic, "^[A-Za-z]+\\s+x?\\s?[A-Za-z]+")
  ) %>% 
  filter(!is.na(species_sci)) %>% # Keep only rows with a valid species_sci
  st_as_sf(
    wkt = "the_geom", # column with WKT strings
    crs = 4326 # EPSG:4326 (lon/lat in WGS84) adjust to local projection for final sample
  ) %>%
  st_filter(denver_county) %>%  # filtering out trees that are not with Denver County
  filter(x_long < -104.7302) %>% 
  left_join(read_csv2("../dc_trees/dc_tree_species.csv"), 
            by = "species_sci") %>% # loading list of dc tree species
  mutate(dc_tree = case_when(
    is.na(dc_tree) ~ FALSE,  # changing na to false
    .default = dc_tree)) %>% 
  left_join(read_csv2("../disease_groups/usda_disease_categories"), 
            by = "disease_pest_def") %>% 
  left_join(read_csv2("../tree_leaf_taxanomy/tree_leaf_taxanomy.csv"), 
            by = "species_sci") %>% 
  mutate(diameter = case_when(diameter == "12 to 1" ~ "12 to 18",
                              diameter == "18 to 2" ~ "18 to 24",
                              diameter == "24 to 3" ~ "24 to 30",
                              diameter == "30 to 3" ~ "30 to 36",
                              diameter == "36 to 4" ~ "36 to 42",
                              diameter == "42 to 4" ~ "42 to 48",
                              .default = diameter)) %>% 
  filter(diameter != "N/A",
         diameter != "0",      
         !is.na(diameter)) %>% 
  filter(!str_detect(species_sci, regex("species", ignore_case = TRUE))) %>% 
  filter(species_sci %in% c("Acer saccharinum", "Ulmus pumila", "Ulmus americana",
                            "Gleditsia triacanthos", "Fraxinus pennsylvanica", 
                            "Populus sargentii", "Pinus nigra",
                            "Populus deltoides", "Catalpa speciosa", 
                            "Celtis occidentalis", "Tilia americana",
                            "Pinus ponderosa", "Quercus rubra", "Ailanthus altissima", 
                            "Tilia cordata"),
         diameter %in% c("24 to 30", "30 to 36", 
                         "36 to 42", "42 to 48", "48 +")) -> denver_trees_subsample

suppressMessages(
  map_dfr(denver_trees_subsample$the_geom, land_cover, .id = "tree_id") %>% 
    mutate(tree_id = as.double(tree_id)))-> lc_data

denver_trees_subsample %>% 
  mutate(tree_id = as.double(row_number())) %>% 
  left_join(lc_data, by = "tree_id") -> denver_trees_subsample

source("../../R/updated_landcover_code.R")

denver_trees_subsample %>%
  rowwise() %>%
  mutate(total = sum(c_across(starts_with("#")), na.rm = TRUE)) %>%
  ungroup() %>% 
  filter(total < 266900) -> mos_sample

suppressMessages(
  map_dfr(mos_sample$the_geom, land_cover_mosaic, .id = "tree_id_2") %>% 
    mutate(tree_id_2 = as.double(tree_id_2))
)-> lc_data

suppressMessages(
  map_dfr(mos_sample$the_geom, land_cover, .id = "tree_id_2") %>% 
    mutate(tree_id_2 = as.double(tree_id_2)))-> lc_data_t

bind_rows(lc_data, lc_data_t) %>% 
  group_by(tree_id_2) %>% 
  summarise(across(starts_with("#"), ~sum(., na.rm = TRUE))) %>%
  rowwise() %>%
  mutate(total = sum(c_across(starts_with("#")), na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(across(starts_with("#_of_pix_"), 
                list(pct = ~. / total), 
                .names = "%_of_{sub('#_of_pix_', '', .col)}")) -> lc_data_t

mos_sample %>% 
  mutate(tree_id_2 = as.double(row_number())) %>% 
  left_join(lc_data_t, by = "tree_id_2") %>% 
  select(-ends_with(".x"), 
         -tree_id_2) %>%
  rename_with(~str_remove(., "\\.y$"), ends_with(".y")) -> mos_sample

denver_trees_subsample %>%
  rowwise() %>%
  mutate(total = sum(c_across(starts_with("#")), na.rm = TRUE)) %>%
  ungroup() %>%
  filter(total < 266900) %>% 
  left_join(mos_sample %>% 
              st_drop_geometry(), by = "tree_id", suffix = c("", "_new")) %>% 
  mutate(across(
    .cols = ends_with("_new"),
    .fns = ~ coalesce(.x, get(sub("_new$", "", cur_column()))),
    .names = "{sub('_new$', '', .col)}"
  )) %>% 
  select(-ends_with("_new"))%>%
  bind_rows(
    denver_trees_subsample %>%
      rowwise() %>%
      mutate(total = sum(c_across(starts_with("#")), na.rm = TRUE)) %>%
      ungroup() %>%
      filter(total >= 266900)
  ) %>%
  select(-starts_with("#"),
         -total) %>% 
  st_write("denver_trees_subsample.geojson", 
           driver = "GeoJSON")
