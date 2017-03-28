#' Melt RNA-Seq count data to long format and log10 transform
#'
#' @author Michael Steinbaugh
#'
#' @import dplyr
#' @import reshape2
#' @import tibble
#'
#' @param bcbio bcbio list object
#' @param counts Counts matrix
#'
#' @return log10 melted data frame
#' @export
melt_log10 <- function(bcbio, counts) {
    metadata <- import_metadata(bcbio) %>%
        dplyr::select_(.dots = unique(c(
            "samplename",
            "description",
            bcbio$intgroup
        )))

    counts <- counts %>%
        as.data.frame %>%
        tibble::rownames_to_column(.) %>%
        reshape2::melt(., id = 1) %>%
        set_names(c("ensembl_gene_id",
                    # Change to `description` instead?
                    "samplename",
                    "counts"))

    # Filter zero counts
    counts <- counts[counts$counts > 0, ]

    # log10 transform
    counts$counts <- log(counts$counts)

    # Join metadata
    counts$samplename <- as.character(counts$samplename)
    counts <- dplyr::left_join(counts, metadata, by = "samplename")

    #str(counts)
    return(counts)
}
