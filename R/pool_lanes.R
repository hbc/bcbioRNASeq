#' Lane-split FASTQ aggregation functions
#'
#' @rdname pool_lanes
#' @author Michael Steinbaugh



#' @rdname pool_lanes
#' @description Pool lane split counts.
#'
#' @param raw_counts Raw counts matrix.
#' @param lane_grep Grep string to match lane identifiers in file name.
#'
#' @return Pooled raw counts matrix.
#' @export
pool_lane_split_counts <- function(raw_counts,
                                   lane_grep = "_L\\d+$") {
    raw_counts <- as.matrix(raw_counts)
    # Obtain the unique pooled sample names
    if (!all(grepl(lane_grep, colnames(raw_counts)))) {
        stop("samples don't appear to be lane split")
    }
    stem <- gsub(lane_grep, "", colnames(raw_counts)) %>% unique %>% sort
    # Perform `rowSums()` on the matching columns per sample
    lapply(seq_along(stem), function(a) {
        raw_counts %>%
            .[, grepl(paste0("^", stem[a], lane_grep), colnames(.))] %>%
            rowSums
    }) %>%
        set_names(stem) %>%
        do.call(cbind, .) %>%
        round
}



#' @rdname pool_lanes
#' @description Pool \linkS4class{DESeqDataSet} counts for technical replicates.
#'   Extract lane split technical replicate counts from a
#'   \linkS4class{DESeqDataSet}, perform [rowSums()], and create a new
#'   \linkS4class{DESeqDataSet} using [DESeq2::DESeqDataSetFromMatrix()].
#'
#' @param run bcbio-nextgen run.
#' @param dds \linkS4class{DESeqDataSet}.
#' @param save Whether to save and export counts.
#'
#' @return \linkS4class{DESeqDataSet} using pooled technical replicates.
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
        normalized_counts <- counts(pooled_dds, normalized = TRUE)
        assign("pooled_normalized_counts",
               normalized_counts,
               envir = parent.frame())
        write_pooled_counts(normalized_counts)

        # raw_counts
        raw_counts <- counts(pooled_dds, normalized = FALSE)
        assign("pooled_raw_counts",
               raw_counts,
               envir = parent.frame())
        write_pooled_counts(raw_counts)
    }

    pooled_dds
}
