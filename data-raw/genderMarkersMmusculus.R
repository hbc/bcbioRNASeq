library(devtools)
load_all()
genderMarkersMmusculus <-
    file.path("data-raw", "genderMarkersMmusculus.csv") %>%
    read_csv
use_data(genderMarkersMmusculus, compress = "xz", overwrite = TRUE)
