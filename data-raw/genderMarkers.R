library(devtools)
library(magrittr)
library(readxl)
load_all()
file <- file.path("data-raw", "genderMarkers.xlsx")
sheets <- excel_sheets(file)
genderMarkers <- lapply(seq_along(sheets), function(a) {
    read_excel(file, sheet = sheets[[a]])
}) %>%
    set_names(sheets)
use_data(genderMarkers, compress = "xz", overwrite = TRUE)
