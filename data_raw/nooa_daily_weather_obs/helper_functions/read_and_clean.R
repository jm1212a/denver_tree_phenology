read_and_clean <- function(file_path) {
  df <- read.fwf(file_path,
                 widths = c(11, 4, 2, 4, rep(8, 31)), 
                 col.names = c("station_id", "year", "month", "element", 
                               paste0("day_", 1:31)),
                 stringsAsFactors = FALSE) %>% 
    mutate(across(starts_with("day_"),
                  ~ as.numeric(str_extract(.x, "-?\\d+"))))
  
  if (nrow(df) == 0) {
    message("Skipping empty file: ", basename(file_path))
    return(NULL)
  }
  if (!all(c("TMAX", "TMIN") %in% df$element)) {
    message("Skipping file missing TMAX or TMIN: ", basename(file_path))
    return(NULL)
  }
  
  df %>%
    mutate(
      year = as.numeric(year),
      month = as.numeric(month)
    ) %>%
    filter(element %in% c("PRCP", "TMAX", "TMIN", "SNWD", "SNOW")) %>%
    pivot_longer(cols = starts_with("day_"), names_to = "day_col", 
                 values_to = "value") %>%
    mutate(
      day = as.numeric(str_extract(day_col, "\\d+")),
      date = ymd(paste(year, month, day, sep = "-")),
      value = na_if(value, 9999),
      value = na_if(value, -9999)
    ) %>%
    filter(!is.na(date)) %>%
    select(-year, -month, -day_col, -day) %>% 
    pivot_wider( names_from = element, values_from = value) %>%
    mutate(TMAX = TMAX / 10,
           TMIN = TMIN / 10)
}
