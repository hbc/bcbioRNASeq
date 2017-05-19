#' Melt RNA-seq count data to long format and log10 transform
#'
#' @keywords internal
#' @author Michael Steinbaugh
#'
#' @param run [bcbioRnaDataSet].
#' @param counts Counts matrix.
#'
#' @return log10 melted data frame.
#' @export
melt_log10 <- function(run, counts) {
    import_tidy_verbs()
    metadata <- run$metadata
    if (!identical(colnames(counts), rownames(metadata))) {
        stop("Sample descriptions in counts do not match metadata.")
    }
    counts %>%
        as.data.frame %>%
        rownames_to_column %>%
        melt(id = 1) %>%
        as_tibble %>%
        set_names(c("ensembl_gene_id",  # rownames
                    "description",  # colnames
                    "counts")) %>%
        # Filter zero counts
        filter(.data$counts > 0) %>%
        # log10 transform
        mutate(counts = log10(.data$counts),
               description = as.character(.data$description)) %>%
        # Join metadata
        left_join(metadata, by = "description")
}



#' Lane-split FASTQ aggregation functions
#'
#' @rdname pool_lanes
#' @author Michael Steinbaugh
#'
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
    # Perform [rowSums()] on the matching columns per sample
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
#' @description Pool [DESeqDataSet] counts for technical replicates. Extract
#'   lane split technical replicate counts from a [DESeqDataSet], perform
#'   [rowSums()], and create a new [DESeqDataSet] using
#'   [DESeqDataSetFromMatrix()].
#'
#' @param run [bcbioRnaDataSet].
#' @param dds [DESeqDataSet].
#' @param save Whether to save and export counts.
#'
#' @return [DESeqDataSet] using pooled technical replicates.
#' @export
pool_lane_split_dds <- function(
    run,
    dds,
    save = TRUE) {
    check_run(run)
    check_dds(dds)

    # Get the internal parameters from [DESeqDataSet]
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

        # Normalized counts
        normalized_counts <- counts(pooled_dds, normalized = TRUE)
        assign("pooled_normalized_counts",
               normalized_counts,
               envir = parent.frame())
        write_pooled_counts(normalized_counts)

        # Raw counts
        raw_counts <- counts(pooled_dds, normalized = FALSE)
        assign("pooled_raw_counts",
               raw_counts,
               envir = parent.frame())
        write_pooled_counts(raw_counts)
    }

    pooled_dds
}



#' Trimmed mean of M-values (TMM) normalization
#'
#' TMM normalization is recommended for RNA-seq data generally when the majority
#' of genes are not differentially expressed. We use this as a quality control
#' tool when plotting counts per gene.
#'
#' @author Michael Steinbaugh
#'
#' @param raw_counts Raw counts matrix.
#'
#' @export
tmm_normalize <- function(raw_counts) {
    raw_counts %>%
        DGEList %>%
        calcNormFactors %>%
        cpm(normalized.lib.sizes = TRUE)
}



#' Transcripts per million
#'
#' Save TPM values from a tximport counts object.
#'
#' @author Michael Steinbaugh
#'
#' @param txi [tximport] list, containing counts in the `abundance`
#'   slot.
#'
#' @return TPM (transcripts per million) matrix.
#' @export
tpm <- function(txi) {
    if (!identical(
        names(txi),
        c("abundance",
          "counts",
          "length",
          "countsFromAbundance")
    )) {
        stop("tximport list is required")
    }
    txi$abundance
}



#' Write gzipped counts
#'
#' @author Michael Steinbaugh
#'
#' @param ... Count matrices, passed in as dots.
#' @param dir Output directory.
#'
#' @export
write_counts <- function(..., dir = file.path("results", "counts")) {
    names <- as.character(substitute(list(...)))[-1L]
    dots <- list(...)

    # Create the counts output directory
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)

    message(paste("Writing", toString(names), "to", dir))

    # Iterate across the dots and write CSVs
    lapply(seq_along(names), function(a) {
        name <- names[a]
        counts <- dots[a]
        counts %>%
            as.data.frame %>%
            rownames_to_column("ensembl_gene_id") %>%
            as_tibble %>%
            write_csv(path = file.path(dir, paste0(name, ".csv.gz")))
    }) %>% invisible
}
