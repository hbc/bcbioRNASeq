#' Pool DESeqDataSet counts for technical replicates
#'
#' Extract lane split technical replicate counts from a \code{DESeqDataSet}
#' object, perform \code{rowSums}, and create a new \code{DESeqDataSet} using
#' \code{DESeqDataSetFromMatrix()}
#'
#' @author Michael Steinbaugh
#'
#' @param run \code{bcbio-nextgen} run object
#' @param dds DESeqDataSet object
#' @param save Whether to save and export counts
#'
#' @return DESeqDataSet object using pooled technical replicates
#' @export
pool_lane_split_dds <- function(
    run,
    dds,
    save = TRUE) {
    check_run(run)
    check_dds(dds)
    name <- deparse(substitute(dds))

    # Get the internal parameters from DESeqDataSet
    raw_counts <- counts(dds, normalized = FALSE)
    design <- DESeq2::design(dds)

    # Pool the lane split technical replicates
    pooled_counts <- pool_lane_split_counts(raw_counts)

    # Obtain metadata
    metadata <- import_metadata(run, pool = TRUE)
    if (!identical(colnames(pooled_counts), rownames(metadata))) {
        stop("count colnames don't match metadata rownames")
    }

    # Re-generate DESeqDataSet using the pooled counts matrix
    dds_pooled <- DESeqDataSetFromMatrix(
        pooled_counts,
        colData = metadata,
        design = design
    ) %>% DESeq

    if (isTRUE(save)) {
        # normalized_counts
        normalized_counts <- counts(dds_pooled, normalized = TRUE)
        assign(paste(name, "pooled_normalized_counts", sep = "_"),
               normalized_counts,
               envir = parent.frame())
        normalized_counts %>%
            rownames_to_column %>%
            write_csv(path = file.path(
                "results",
                paste(name, "pooled_normalized_counts.csv.gz", sep = "_")))

        # raw_counts
        raw_counts <- counts(dds_pooled, normalized = FALSE)
        assign(paste(name, "pooled_raw_counts", sep = "_"),
               raw_counts,
               envir = parent.frame())
        raw_counts %>%
            rownames_to_column %>%
            write_csv(path = file.path(
                "results",
                paste(name, "pooled_raw_counts.csv.gz", sep = "_")))
    }

    return(dds_pooled)
}
