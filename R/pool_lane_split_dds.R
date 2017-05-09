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

    # Get the internal parameters from DESeqDataSet
    raw_counts <- counts(dds, normalized = FALSE)
    design <- design(dds)

    # Pool the lane split technical replicates
    pooled_counts <- pool_lane_split_counts(raw_counts)

    # Reload metadata without lanes
    if (is.null(run$metadata_file)) {
        stop("Custom metadata file required for count pooling")
    }
    metadata <- read_metadata(run$metadata_file) %>%
        .[order(.$file_name), ]
    if (!identical(colnames(pooled_counts), metadata$file_name)) {
        stop("Pooled count names don't match original file names")
    }

    # Replace file names in pooled count matrix with description
    colnames(pooled_counts) <- metadata$description

    # Ensure count matrix colnames matches metadata rownames
    if (!identical(colnames(pooled_counts), rownames(metadata))) {
        stop("Description mismatch")
    }

    pooled_dds <- DESeqDataSetFromMatrix(
        countData = pooled_counts,
        colData = metadata,
        design = design
    ) %>% DESeq

    if (isTRUE(save)) {
        # Internal rework of `write_counts()` that doesn't use dots
        write_pooled_counts <- function(counts, dir = "results/counts") {
            name <- deparse(substitute(counts))
            counts %>%
                as.data.frame %>%
                rownames_to_column("ensembl_gene_id") %>%
                as_tibble %>%
                write_csv(path = file.path(dir, paste0(
                    paste("pooled", name, sep = "_"), ".csv.gz")))
        }

        # normalized_counts
        normalized_counts <- counts(dds_pooled, normalized = TRUE)
        assign("pooled_normalized_counts",
               normalized_counts,
               envir = parent.frame())
        write_pooled_counts(normalized_counts)

        # raw_counts
        raw_counts <- counts(dds_pooled, normalized = FALSE)
        assign("pooled_raw_counts",
               raw_counts,
               envir = parent.frame())
        write_pooled_counts(raw_counts)
    }

    return(pooled_dds)
}
