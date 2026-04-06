library(tidyverse)

fix_hex_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  
  # Recover 24eN from 2.40E+N pattern (handles both 2.4e+N and 2.40e+N)
  sci_idx <- grepl("^2\\.40?e\\+0?([0-9]+)$", x, ignore.case = TRUE)
  x[sci_idx] <- paste0("24e", sub("^2\\.40?e\\+0?([0-9]+)$", "\\1",
                                  x[sci_idx], ignore.case = TRUE))
  
  # Null out anything else still in scientific notation
  x[grepl("^[0-9.]+e[+-][0-9]+$", x, ignore.case = TRUE) & !sci_idx] <- NA_character_
  x[tolower(x) %in% c("0", "na", "")] <- NA_character_
  
  tolower(x)
}