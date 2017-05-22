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
#' Pool lane-split counts matrix or [DESeqDataSet] using [rowSums()].
#'
#' @rdname pool_lanes
#' @author Michael Steinbaugh
#'
#' @param counts Raw counts matrix.
#' @param pattern Grep pattern to match lane identifiers in file name.
#'
#' @return Object of same class, with pooled technical replicates.
#' @export
pool_counts <- function(counts, pattern = "_L\\d+") {
    # Ensure that counts are input as matrix
    counts <- as.matrix(counts)
    # Obtain the unique pooled sample names
    if (!all(grepl(pattern, colnames(counts)))) {
        stop("Samples don't match lane grep pattern")
    }
    stem <- str_replace(colnames(counts), pattern, "") %>% unique %>% sort
    # Perform [rowSums()] on the matching columns per sample
    lapply(seq_along(stem), function(a) {
        counts %>%
            .[, grepl(paste0("^", stem[a], pattern), colnames(.))] %>%
            rowSums
    }) %>%
        set_names(stem) %>%
        do.call(cbind, .) %>%
        # [round()] here otherwise [DESeq()] will error out
        round
}



#' @rdname pool_lanes
#' @param dds [DESeqDataSet].
#' @param save Whether to assign into environment and write CSV files.
#' @export
#' @seealso A new [DESeqDataSet] is returned from [DESeqDataSetFromMatrix()].
pool_dds <- function(dds, pattern = "_L\\d+", save = FALSE) {
    check_dds(dds)
    import_tidy_verbs()
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
