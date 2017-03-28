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
    metadata <- import_metadata(bcbio)

    counts <- tmm_normalized_counts
    counts <- as.data.frame(counts)
    counts <- tibble::rownames_to_column(counts)
    counts <- reshape2::melt(counts, id = 1)
    counts <- set_names(counts, c("ensembl_gene_id",
                                  # Change to `description` instead?
                                  "samplename",
                                  "counts"))
    # Filter zero counts
    counts <- counts[counts$counts > 0, ]
    # log10 transform
    counts$counts <- log(counts$counts)
    # Convert samplename from factor to character
    counts$samplename <- as.character(counts$samplename)
    # Join metadata
    counts <- dplyr::left_join(counts, metadata, by = "samplename")

    #str(counts)
    return(counts)
}
