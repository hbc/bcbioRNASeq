library(devtools)
library(readxl)
library(tidyverse)
load_all()
file <- "data-raw/genderMarkers.xlsx"
sheets <- excel_sheets(file)
genderMarkers <- lapply(seq_along(sheets), function(a) {
    read_excel(file, sheet = sheets[[a]])
})
names(genderMarkers) <- camel(sheets)
use_data(genderMarkers, compress = "xz", overwrite = TRUE)
