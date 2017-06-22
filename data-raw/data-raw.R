library(basejump)
gender_markers_mmusculus <-
    file.path("data-raw", "gender_markers_mmusculus.csv") %>%
    read_csv
save_data(gender_markers_mmusculus)
load_all()
