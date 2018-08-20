library(devtools)
library(readxl)
library(tidyverse)

file <- "data-raw/gender_markers.xlsx"
sheets <- excel_sheets(file)
gender_markers <- lapply(seq_along(sheets), function(a) {
    read_excel(file, sheet = sheets[[a]]) %>%
        filter(include == TRUE) %>%
        select(-include)
})
names(gender_markers) <- camel(sheets)
use_data(gender_markers, compress = "xz", overwrite = TRUE)
