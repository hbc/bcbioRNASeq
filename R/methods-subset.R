#' Bracket-Based Subsetting
#'
#' Extract genes by row and samples by column from a `bcbioRNASeq` object. The
#' internal `DESeqDataSet` and count transformations are rescaled automatically.
#' DESeq2 transformations can be disabled on large subset operations by setting
#' `transform = FALSE`.
#'
#' @name subset
#' @family Data Functions
#' @author Lorena Pantano, Michael Steinbaugh
#'
#' @importFrom DESeq2 DESeq DESeqDataSetFromTximport rlog
#'   varianceStabilizingTransformation
#'
#' @inheritParams base::`[`
#' @param ... Additional arguments.
#'
#' @return `bcbioRNASeq`.
#'
#' @seealso `help("[", "base")`.
#'
#' @examples
#' genes <- rownames(bcb_small)[1:50]
#' head(genes)
#' samples <- colnames(bcb_small)[1:2]
#' head(samples)
#'
#' # Subset by sample name
#' bcb_small[, samples]
#'
#' # Subset by gene list
#' bcb_small[genes, ]
#'
#' # Subset by both genes and samples
#' subset <- bcb_small[genes, samples]
#' print(subset)
#' assayNames(subset)
#'
#' # Skip DESeq2 variance stabilizing transformations
#' subset <- bcb_small[genes, samples, transform = FALSE]
#' print(subset)
#' assayNames(subset)
NULL



# Constructors =================================================================
#' @importFrom DESeq2 DESeq estimateSizeFactors rlog
#'   varianceStabilizingTransformation
#' @importFrom tibble column_to_rownames rownames_to_column
.subset.bcbioRNASeq <- function(x, i, j, ..., drop) {  # nolint
    validObject(x)

    dots <- list(...)
    if (!identical(dots[["transform"]], FALSE)) {
        transform <- TRUE
    }

    # Genes (rows)
    if (missing(i)) {
        i <- 1L:nrow(x)
    }
    # Require at least 100 genes
    assert_all_are_in_left_open_range(length(i), lower = 99L)

    # Samples (columns)
    if (missing(j)) {
        j <- 1L:ncol(x)
    }
    # Require at least 2 samples
    assert_all_are_in_left_open_range(length(j), lower = 1L)


    # Early return if dimensions are unmodified
    if (identical(dim(x), c(length(i), length(j)))) {
        return(x)
    }

    # Regenerate RangedSummarizedExperiment
    rse <- as(x, "RangedSummarizedExperiment")
    rse <- rse[i, j, drop = drop]

    # Row and column data ======================================================
    rowRanges <- rowRanges(rse)
    colData <- colData(rse)

    # Assays ===================================================================
    assays <- assays(rse)

    if (isTRUE(transform)) {
        inform(paste(
            "Calculating variance stabilizations using DESeq2",
            packageVersion("DESeq2")
        ))
        txi <- .regenerateTximportList(rse)
        dds <- DESeqDataSetFromTximport(
            txi = txi,
            colData = colData,
            # Use an empty design formula
            design = ~ 1  # nolint
        )
        # Suppress warning about empty design formula
        dds <- suppressWarnings(DESeq(dds))
        validObject(dds)
        # FIXME This isn't using the correct dimensions
        inform("Applying rlog transformation")
        assays[["rlog"]] <- assay(rlog(dds))
        inform("Applying variance stabilizing transformation")
        assays[["vst"]] <- assay(varianceStabilizingTransformation(dds))
    } else {
        inform("Skipping DESeq2 transformations")
        # Ensure existing transformations get dropped
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
        .[colnames(rse), , drop = FALSE] %>%
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
        j = "ANY",
        drop = "ANY"  # don't use logical here
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
