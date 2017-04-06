#' Pool lane split counts
#'
#' @author Michael Steinbaugh
#'
#' @param counts Counts matrix
#' @param lane_grep Grep string to match lane identifiers in file name
#'
#' @return Pooled counts matrix
#' @export
pool_lane_split_counts <- function(counts,
                                   lane_grep = "_L\\d+$") {
    counts <- as.matrix(counts)
    # Obtain the unique pooled sample names
    if (!all(grepl(lane_grep, colnames(counts)))) {
        stop("samples don't appear to be lane split")
    }
    stem <- gsub(lane_grep, "", colnames(counts)) %>% unique %>% sort
    # Perform `rowSums()` on the matching columns per sample
    pooled_counts <- lapply(seq_along(stem), function(a) {
        counts %>%
            .[, grepl(paste0("^", stem[a], lane_grep), colnames(.))] %>%
            rowSums
    }) %>%
        set_names(stem) %>%
        do.call(cbind, .) %>%
        # Counts must be integers!
        round

    return(pooled_counts)
}

