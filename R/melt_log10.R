#' Melt RNA-Seq count data to long format and log10 transform
#'
#' @author Michael Steinbaugh
#'
#' @param run \code{bcbio-nextgen} run
#' @param counts Counts matrix
#'
#' @return log10 melted data frame
#' @export
melt_log10 <- function(run, counts) {
    check_run(run)
    metadata <- run$metadata
    if (!identical(colnames(counts), rownames(metadata))) {
        stop("Sample descriptions in counts do not match metadata.")
    }

    melted <- counts %>%
        as.data.frame %>%
        rownames_to_column %>%
        melt(id = 1) %>%
        as_tibble %>%
        set_names(c("ensembl_gene_id",  # rownames
                    "description",  # colnames
                    "counts")) %>%
        # Filter zero counts
        .[.$counts > 0, ] %>%
        # log10 transform
        mutate(counts = log10(.data$counts),
               description = as.character(.data$description)) %>%
        # Join metadata
        left_join(metadata, by = "description")

    return(melted)
}
