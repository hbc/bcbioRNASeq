devtools::load_all()
library(tidyverse)
file <- "data-raw/genderMarkers.xlsx"
sheets <- excel_sheets(file)
genderMarkers <- lapply(seq_along(sheets), function(a) {
    read_excel(file, sheet = sheets[[a]])
})
names(genderMarkers) <- camel(sheets)
devtools::use_data(genderMarkers, compress = "xz", overwrite = TRUE)
