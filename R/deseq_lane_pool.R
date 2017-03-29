#' Run DESeq2 with pooled technical replicates
#'
#' @author Michael Steinbaugh
#'
#' @import DESeq2
#'
#' @param bcbio bcbio list object
#' @param dds DESeqDataSet object
#'
#' @return DESeqDataSet object using pooled technical replicates
#' @export
deseq_lane_pool <- function(bcbio, dds) {
    if (class(dds)[1] != "DESeqDataSet") {
        stop("A DESeqDataSet is required.")
    }

    # Get the internal parameters from DESeqDataSet
    counts <- DESeq2::counts(dds)
    design <- DESeq2::design(dds)

    # Obtain the unique pooled sample names
    lane_grep <- "_L\\d+$"
    stem <- gsub(lane_grep, "", colnames(counts)) %>% unique %>% sort

    # Perform `rowSums` on the matching columns per sample
    pooled_counts <- lapply(seq_along(stem), function(a) {
        counts %>%
            .[, grepl(paste0("^", stem[a], lane_grep), colnames(.))] %>%
            rowSums
    }) %>%
        set_names(stem) %>%
        do.call(cbind, .) %>%
        # Counts must be integers!
        round

    # Obtain lane pool metadata
    metadata <- import_metadata(bcbio, lanes = "pooled")

    if (!identical(colnames(pooled_counts), rownames(metadata))) {
        stop("Count column names don't match the metadata row names.")
    }

    # Run DESeq2 using the pooled counts matrix
    dds_pooled <- DESeq2::DESeqDataSetFromMatrix(
        pooled_counts,
        colData = metadata,
        design = design
    ) %>% DESeq2::DESeq(.)

    return(dds_pooled)
}
