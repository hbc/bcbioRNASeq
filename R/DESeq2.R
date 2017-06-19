#' @rdname pool
#' @param dds [DESeqDataSet].
#' @param save Whether to assign into environment and write CSV files.
#' @export
#' @seealso A new [DESeqDataSet] is returned from [DESeqDataSetFromMatrix()].
pool_dds <- function(dds, pattern = "_L\\d+", save = FALSE) {
    check_dds(dds)
    envir <- parent.frame()

    # Get the stored design formula
    design <- design(dds)

    # Pool the lane split technical replicates
    countData <- counts(dds, normalized = FALSE) %>% pool_counts

    # Mutate the colData metadata to pooled samples
    colData <- colData(dds) %>%
        as.data.frame %>%
        mutate(file_name = str_replace(.data$file_name, pattern, ""),
               description = str_replace(.data$description, pattern, ""),
               lane = NULL,
               # [DESeq()]-generated columns
               # replaceable = NULL,
               sizeFactor = NULL) %>%
        distinct %>%
        set_rownames(.$description)

    # Check that the new colData matches the counts matrix
    if (!identical(colData$file_name, colnames(countData))) {
        stop("File name mismatch in colData and countData")
    }

    # Replace file names in pooled count matrix with description
    colnames(countData) <- colData$description
    if (!identical(rownames(colData), colnames(countData))) {
        stop("Description mismatch in colData and countData")
    }

    # Pipe back into DESeq
    pooled_dds <- DESeqDataSetFromMatrix(
        countData = countData,
        colData = colData,
        design = design
    ) %>% DESeq

    # Write CSVs and assign pooled counts to environment, if desired
    if (isTRUE(save)) {
        # Internal rework of [write_counts()] that doesn't use dots
        write_pooled_counts <- function(counts, dir = "results/counts") {
            dir.create(dir, recursive = TRUE, showWarnings = FALSE)
            name <- deparse(substitute(counts))
            counts %>%
                as.data.frame %>%
                rownames_to_column("ensembl_gene_id") %>%
                as_tibble %>%
                write_csv(path = file.path(dir, paste0(
                    paste("pooled", name, sep = "_"), ".csv.gz")))
        }

        # Normalized counts
        normalized_counts <- counts(pooled_dds, normalized = TRUE)
        assign("pooled_normalized_counts", normalized_counts, envir = envir)
        write_pooled_counts(normalized_counts)

        # Raw counts
        raw_counts <- counts(pooled_dds, normalized = FALSE)
        assign("pooled_raw_counts", raw_counts, envir = envir)
        write_pooled_counts(raw_counts)
    }

    pooled_dds
}
