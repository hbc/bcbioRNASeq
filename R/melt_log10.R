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
    metadata <- import_metadata(run) %>%
        select_(.dots = unique(c(
            "samplename",
            "description",
            run$intgroup
        )))

    if (!identical(colnames(counts), rownames(metadata))) {
        stop("Sample descriptions in counts do not match metadata.")
    }

    melted <- counts %>%
        as.data.frame %>%
        rownames_to_column %>%
        melt(id = 1) %>%
        set_names(c("ensembl_gene_id",  # rownames
                    "description",  # colnames
                    "counts"))

    # Filter zero counts
    melted <- melted[melted$counts > 0, ]

    # log10 transform
    melted$counts <- log(melted$counts)

    # Join metadata
    melted$description <- as.character(melted$description)
    melted <- left_join(melted, metadata, by = "description")
    #str(melted)

    return(melted)
}
