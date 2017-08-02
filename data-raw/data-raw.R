library(basejump)
genderMarkersMmusculus <-
    file.path("data-raw", "genderMarkersMmusculus.csv") %>%
    read_csv
saveData(genderMarkersMmusculus)
devtools::load_all()
