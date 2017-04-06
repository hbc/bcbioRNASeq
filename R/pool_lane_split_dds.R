#' Pool DESeqDataSet counts for technical replicates
#'
#' Extract lane split technical replicate counts from a \code{DESeqDataSet}
#' object, perform \code{rowSums}, and create a new \code{DESeqDataSet} using
#' \code{DESeqDataSetFromMatrix()}
#'
#' @author Michael Steinbaugh
#'
#' @import DESeq2
#'
#' @param bcbio bcbio list object
#' @param dds DESeqDataSet object
#' @param save_counts Whether to save and export counts
#'
#' @return DESeqDataSet object using pooled technical replicates
#' @export
pool_lane_split_dds <- function(
    bcbio,
    dds,
    save_counts = TRUE) {
    check_bcbio_object(bcbio)
    if (class(dds)[1] != "DESeqDataSet") {
        stop("A DESeqDataSet is required.")
    }
    name <- deparse(substitute(dds))

    # Get the internal parameters from DESeqDataSet
    raw_counts <- DESeq2::counts(dds, normalized = FALSE)
    design <- DESeq2::design(dds)

    # Pool the lane split technical replicates
    pooled_counts <- pool_lane_split_counts(raw_counts)

    # Obtain metadata
    metadata <- import_metadata(bcbio, pool = TRUE)
    if (!identical(colnames(pooled_counts), rownames(metadata))) {
        stop("Count column names don't match the metadata row names.")
    }

    # Re-generate DESeqDataSet using the pooled counts matrix
    dds_pooled <- DESeq2::DESeqDataSetFromMatrix(
        pooled_counts,
        colData = metadata,
        design = design
    ) %>% DESeq2::DESeq(.)

    if (isTRUE(save_counts)) {
        # normalized_counts
        normalized_counts <- DESeq2::counts(dds_pooled, normalized = TRUE)
        assign(paste(name, "pooled_normalized_counts", sep = "_"),
               normalized_counts,
               envir = parent.frame())
        write.csv(
            normalized_counts,
            file = file.path(
                "results",
                paste(name,
                      "pooled_normalized_counts.csv", sep = "_")))

        # raw_counts
        raw_counts <- DESeq2::counts(dds_pooled, normalized = FALSE)
        assign(paste(name, "pooled_raw_counts", sep = "_"),
               raw_counts,
               envir = parent.frame())
        write.csv(
            raw_counts,
            file = file.path(
                "results",
                paste(name,
                      "pooled_raw_counts.csv", sep = "_")))
    }

    return(dds_pooled)
}
