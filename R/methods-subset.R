#' Bracket-Based Subsetting
#'
#' Extract genes by row and samples by column from a `bcbioRNASeq` object. The
#' internal `DESeqDataSet` and count transformations are rescaled automatically.
#' DESeq2 transformations can be disabled on large subset operations by setting
#' `transform = FALSE`.
#'
#' @name subset
#' @author Lorena Pantano, Michael Steinbaugh
#'
#' @inheritParams base::`[`
#' @param ... Additional arguments.
#'
#' @return `bcbioRNASeq`.
#'
#' @seealso `help("[", "base")`.
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))
#'
#' genes <- rownames(bcb)[1:50]
#' head(genes)
#' samples <- colnames(bcb)[1:2]
#' head(samples)
#'
#' # Subset by sample name
#' bcb[, samples]
#'
#' # Subset by gene list
#' bcb[ensgene, ]
#'
#' # Subset by both genes and samples
#' bcb[ensgene, samples]
#'
#' # Skip DESeq2 transformations
#' bcb[ensgene, samples, transform = FALSE]
NULL



# Constructors =================================================================
#' @importFrom DESeq2 DESeq estimateSizeFactors rlog
#'   varianceStabilizingTransformation
#' @importFrom tibble column_to_rownames rownames_to_column
.subset.bcbioRNASeq <- function(x, i, j, ..., drop = FALSE) {  # nolint
    validObject(x)

    # Genes (rows)
    if (missing(i)) {
        i <- 1L:nrow(x)
    }
    assert_all_are_greater_than(length(i), 1L)

    # Samples (columns)
    if (missing(j)) {
        j <- 1L:ncol(x)
    }
    assert_all_are_greater_than(length(j), 1L)

    # Early return if dimensions are unmodified
    if (identical(dim(x), c(length(i), length(j)))) {
        return(x)
    }

    dots <- list(...)
    if (!identical(dots[["transform"]], FALSE)) {
        transform <- TRUE
    }

    # Regenerate RangedSummarizedExperiment
    rse <- as(x, "RangedSummarizedExperiment")
    rse <- rse[i, j, drop = drop]

    genes <- rownames(rse)
    samples <- colnames(rse)

    # Row ranges ===============================================================
    rowRanges <- rowRanges(rse)

    # Column data ==============================================================
    colData <- colData(rse)
    colData <- sanitizeColData(colData)

    # Assays ===================================================================
    assays <- assays(rse)

    # Update the colData in nested DESeq objects
    message("Updating colData in assays slot")
    colData(assays[["dds"]]) <- colData
    if (is(assays[["rlog"]], "DESeqTransform")) {
        colData(assays[["rlog"]]) <- colData
    }
    if (is(assays[["vst"]], "DESeqTransform")) {
        colData(assays[["vst"]]) <- colData
    }

    # Update DESeqTransform objects
    if (identical(transform, TRUE)) {
        inform("Updating internal DESeqDataSet")
        # DESeq2 will warn about empty design formula, if set
        assays[["dds"]] <- suppressWarnings(DESeq(assays[["dds"]]))
        inform("Performing rlog transformation")
        assays[["rlog"]] <- rlog(assays[["dds"]])
        inform("Performing variance stabilizing transformation")
        assays[["vst"]] <- varianceStabilizingTransformation(assays[["dds"]])
    } else if (identical(transform, FALSE)) {
        # Skip normalization option, for large datasets.
        # Not generally recommended.
        inform("Updating size factors for internal DESeqDataset")
        assays[["dds"]] <- estimateSizeFactors(assays[["dds"]])
        inform("Skipping DESeq2 transformations")
        assays[["vst"]] <- NULL
        assays[["rlog"]] <- NULL
    }

    # Drop any NULL items in assays list
    assays <- Filter(Negate(is.null), assays)

    # Metadata =================================================================
    metadata <- metadata(x)
    metadata[["subset"]] <- TRUE
    # Update version, if necessary
    if (!identical(metadata[["version"]], packageVersion)) {
        metadata[["originalVersion"]] <- metadata[["version"]]
        metadata[["version"]] <- packageVersion
    }
    # Metrics
    metadata[["metrics"]] <- metadata[["metrics"]] %>%
        .[samples, , drop = FALSE] %>%
        rownames_to_column() %>%
        mutate_if(is.character, as.factor) %>%
        mutate_if(is.factor, droplevels) %>%
        column_to_rownames()

    # Return ===================================================================
    rse <- SummarizedExperiment(
        assays = assays,
        rowRanges = rowRanges,
        colData = colData,
        metadata = metadata
    )
    assert_is_all_of(rse, "RangedSummarizedExperiment")
    new("bcbioRNASeq", rse)
}



# Methods ======================================================================
#' @rdname subset
#' @export
setMethod(
    "[",
    signature(
        x = "bcbioRNASeq",
        i = "ANY",
        j = "ANY"
    ),
    function(x, i, j, ..., drop = FALSE) {
        .subset.bcbioRNASeq(
            x = x,
            i = i,
            j = j,
            ...,
            drop = drop
        )
    }
)
