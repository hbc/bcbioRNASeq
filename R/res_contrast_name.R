res_contrast_name <- function(res) {
    mcols(res)[2, 2] %>%
        gsub("^.*:\\s", "", .)
}
