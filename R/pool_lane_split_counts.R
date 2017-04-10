#' Pool lane split counts
#'
#' @author Michael Steinbaugh
#'
#' @param raw_counts Raw counts matrix
#' @param lane_grep Grep string to match lane identifiers in file name
#'
#' @return Pooled raw counts matrix
#' @export
pool_lane_split_counts <- function(raw_counts,
                                   lane_grep = "_L\\d+$") {
    raw_counts <- as.matrix(raw_counts)
    # Obtain the unique pooled sample names
    if (!all(grepl(lane_grep, colnames(raw_counts)))) {
        stop("samples don't appear to be lane split")
    }
    stem <- gsub(lane_grep, "", colnames(raw_counts)) %>% unique %>% sort
    # Perform `rowSums()` on the matching columns per sample
    pooled_counts <- lapply(seq_along(stem), function(a) {
        raw_counts %>%
            .[, grepl(paste0("^", stem[a], lane_grep), colnames(.))] %>%
            rowSums
    }) %>%
        set_names(stem) %>%
        do.call(cbind, .) %>%
        round

    return(pooled_counts)
}

