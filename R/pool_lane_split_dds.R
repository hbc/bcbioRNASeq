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
#' @param assign_counts Whether to save and export counts
#'
#' @return DESeqDataSet object using pooled technical replicates
#' @export
pool_lane_split_dds <- function(
    bcbio,
    dds,
    assign_counts = TRUE) {
    check_bcbio_object(bcbio)
    if (class(dds)[1] != "DESeqDataSet") {
        stop("A DESeqDataSet is required.")
    }

    name <- deparse(substitute(dds))

    # Get the internal parameters from DESeqDataSet
    raw_counts <- DESeq2::counts(dds, normalized = FALSE)
    design <- DESeq2::design(dds)

    # Obtain the unique pooled sample names
    lane_grep <- "_L\\d+$"
    if (!all(grepl(lane_grep, colnames(raw_counts)))) {
        stop("Samples aren't lane split, and don't need to be pooled.")
    }
    stem <- gsub(lane_grep, "", colnames(raw_counts)) %>%
        unique %>% sort

    # Perform `rowSums` on the matching columns per sample
    pooled_counts <- lapply(seq_along(stem), function(a) {
        raw_counts %>%
            .[, grepl(paste0("^", stem[a], lane_grep), colnames(.))] %>%
            rowSums
    }) %>%
        set_names(stem) %>%
        do.call(cbind, .) %>%
        # Counts must be integers!
        round
    rm(raw_counts)

    # Obtain lane pool metadata
    metadata <- import_metadata(bcbio, pool = TRUE)
    if (!identical(colnames(pooled_counts), rownames(metadata))) {
        stop("Count column names don't match the metadata row names.")
    }

    # Run DESeq2 using the pooled counts matrix
    dds_pooled <- DESeq2::DESeqDataSetFromMatrix(
        pooled_counts,
        colData = metadata,
        design = design
    ) %>% DESeq2::DESeq(.)

    if (isTRUE(assign_counts)) {
        # normalized_counts
        normalized_counts <- DESeq2::counts(dds_pooled, normalized = TRUE)
        assign(paste(name, "pooled_normalized_counts", sep = "_"),
               normalized_counts,
               envir = parent.frame())
        # raw_counts
        raw_counts <- DESeq2::counts(dds_pooled, normalized = FALSE)
        assign(paste(name, "pooled_raw_counts", sep = "_"),
               raw_counts,
               envir = parent.frame())
    }

    return(dds_pooled)
}
