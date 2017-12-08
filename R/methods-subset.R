#' Bracket-Based Subsetting
#'
#' Extract genes by row and samples by column from a [bcbioRNASeq] object. The
#' internal [DESeqDataSet] and count transformations are rescaled automatically.
#'
#' @rdname subset
#' @name subset
#'
#' @author Lorena Pantano, Michael Steinbaugh
#'
#' @inheritParams base::`[`
#' @param ... Additional arguments.
#'
#' @return [bcbioRNASeq].
#'
#' @seealso `help("[", "base")`.
#'
#' @examples
#' load(system.file(
#'     file.path("inst", "extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
#'
#' ensgene <- rownames(bcb)[1:50]
#' head(ensgene)
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
#' # Skip normalization
#' bcb[ensgene, samples, skipNorm = TRUE]
NULL



# Constructors ====
#' @importFrom DESeq2 DESeq estimateSizeFactors rlog
#'   varianceStabilizingTransformation
#' @importFrom dplyr mutate_if
#' @importFrom S4Vectors metadata SimpleList
#' @importFrom tibble column_to_rownames rownames_to_column
.subset <- function(x, i, j, ..., drop = FALSE) {
    # Genes (rows)
    if (missing(i)) {
        i <- 1:nrow(x)
    }
    if (length(i) < 2) {
        stop("At least 2 genes are required", call. = FALSE)
    }

    # Samples (columns)
    if (missing(j)) {
        j <- 1:ncol(x)
    }
    if (length(j) < 2) {
        stop("At least 2 samples are required", call. = FALSE)
    }

    # Early return if dimensions are unmodified
    if (identical(dim(x), c(length(i), length(j)))) return(x)

    dots <- list(...)
    if (is.null(dots[["transformationLimit"]])) {
        transformationLimit <- 50
    } else {
        transformationLimit <- dots[["transformationLimit"]]
    }
    if (is.null(dots[["skipNorm"]])) {
        skipNorm <- FALSE
    } else {
        skipNorm <- dots[["skipNorm"]]
    }

    # Regenerate and subset SummarizedExperiment
    se <- as(x, "SummarizedExperiment")
    se <- se[i, j, drop = drop]

    genes <- rownames(se)
    samples <- colnames(se)

    # Row data =================================================================
    rowData <- rowData(se)
    if (!is.null(rowData)) {
        rownames(rowData) <- slot(se, "NAMES")
    }

    # Column data ==============================================================
    # Better base R approach here to relevel factors?
    colData <- colData(se) %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        mutate_if(is.character, as.factor) %>%
        mutate_if(is.factor, droplevels) %>%
        column_to_rownames() %>%
        as("DataFrame")

    # bcbio ====================================================================
    # tximport
    txi <- bcbio(x, "tximport")
    txi[["abundance"]] <- txi[["abundance"]][genes, samples]
    txi[["counts"]] <- txi[["counts"]][genes, samples]
    txi[["length"]] <- txi[["length"]][genes, samples]

    # DESeqDataSet
    message("Updating internal DESeqDataSet")
    dds <- bcbio(x, "DESeqDataSet")
    dds <- dds[genes, samples]
    colData(dds) <- colData
    # Skip normalization option, for large datasets
    if (isTRUE(skipNorm) | nrow(colData) > transformationLimit) {
        message("Skipping re-normalization, just selecting samples and genes")
        dds <- estimateSizeFactors(dds)
        vst <- NULL
        rlog <- NULL
    } else {
        # DESeq2 will warn about empty design formula, if set
        dds <- suppressWarnings(DESeq(dds))
        message("Performing rlog transformation")
        rlog <- rlog(dds)
        message("Performing variance stabilizing transformation")
        vst <- varianceStabilizingTransformation(dds)
    }

    # featureCounts
    featureCounts <- bcbio(x, "featureCounts")
    if (is.matrix(featureCounts)) {
        featureCounts <- featureCounts[genes, samples, drop = FALSE]
    } else {
        featureCounts <- NULL
    }

    bcbio <- SimpleList(
        tximport = txi,
        DESeqDataSet = dds,
        featureCounts = featureCounts)

    # Assays ===================================================================
    raw <- txi[["counts"]]
    normalized <- counts(dds, normalized = TRUE)
    tpm <- txi[["abundance"]]
    tmm <- tmm(raw)
    assays <- SimpleList(
        raw = raw,
        normalized = normalized,
        tpm = tpm,
        tmm = tmm,
        rlog = rlog,
        vst = vst)
    # Drop `NULL` assay slots, if necessary. We need this step if rlog and vst
    # transformations are skipped above.
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
        .[.[["sampleID"]] %in% samples, , drop = FALSE] %>%
        rownames_to_column() %>%
        mutate_if(is.character, as.factor) %>%
        mutate_if(is.factor, droplevels) %>%
        column_to_rownames()

    # Return ===================================================================
    new("bcbioRNASeq",
        SummarizedExperiment(
            assays = assays,
            rowData = rowData,
            colData = colData,
            metadata = metadata),
        bcbio = bcbio)
}



# Methods ====
#' @rdname subset
#' @export
setMethod(
    "[",
    signature(x = "bcbioRNASeq",
              i = "ANY",
              j = "ANY"),
    function(x, i, j, ..., drop = FALSE) {
        .subset(x, i, j, ..., drop)
    })
