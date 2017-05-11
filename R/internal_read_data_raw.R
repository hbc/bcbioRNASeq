#' Read and save files in data-raw directory
#'
#' @keywords internal
#' @author Michael Steinbaugh
#'
#' @export
read_data_raw <- function() {
    csv <- list.files("data-raw", pattern = "*.csv", full.names = TRUE)
    sapply(seq_along(csv), function(a) {
        name <- basename(csv[a]) %>% file_path_sans_ext
        df <- read_csv(csv[a])
        dir.create("data", showWarnings = FALSE)
        assign(name, df)
        save(list = name, file = file.path("data", paste0(name, ".rda")))
    }) %>% invisible
}
