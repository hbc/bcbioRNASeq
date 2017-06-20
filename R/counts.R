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
melt_log10 <- function(bcb, counts = NULL) {
    metadata <- colData(bcb) %>% as.data.frame
    if (is.null(counts)) {
        counts <- assay(bcb)
    }
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
        filter(.data[["counts"]] > 0) %>%
        # log10 transform
        mutate(counts = log10(.data[["counts"]]),
               description = as.character(.data[["description"]])) %>%
        # Join metadata
        left_join(metadata, by = "description")
}



#' Lane-split FASTQ aggregation functions
#'
#' Pool lane-split counts matrix or [DESeqDataSet] using [rowSums()].
#'
#' @rdname pool
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
