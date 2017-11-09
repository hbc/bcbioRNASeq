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
#' data(bcb)
#' genes <- 1:50
#' samples <- c("group1_1", "group1_2")
#'
#' # Subset by sample name
#' bcb[, samples]
#'
#' # Subset by gene list
#' bcb[genes, ]
#'
#' # Subset by both genes and samples
#' \dontrun{
#' bcb[genes, samples]
#' }
#'
#' # Skip normalization
#' \dontrun{
#' bcb[genes, samples, skipNorm = TRUE]
#' }
NULL



# Constructors ====
#' Create DDS
#'
#' This operation must be placed outside of the S4 method dispatch. Otherwise,
#' the resulting subset object will be ~2X the expected size on disk when
#' saving, for an unknown reason.
#'
#' @noRd
#'
#' @importFrom DESeq2 DESeqDataSetFromTximport
#' @importFrom stats formula
.createDDS <- function(txi, colData) {
    DESeqDataSetFromTximport(
        txi = txi,
        colData = colData,
        design = formula(~1))
}



#' @importFrom DESeq2 DESeqTransform
.subsetCounts <- function(counts, colData) {
    se <- SummarizedExperiment(assays = counts, colData = colData)
    DESeqTransform(se)
}



#' @importFrom DESeq2 DESeq estimateSizeFactors rlog
#'   varianceStabilizingTransformation
#' @importFrom S4Vectors metadata SimpleList
.subset <- function(x, i, j, ..., drop = FALSE) {
    if (missing(i)) {
        i <- 1:nrow(x)
    }
    if (missing(j)) {
        j <- 1:ncol(x)
    }

    dots <- list(...)
    if (is.null(dots[["maxSamples"]])) {
        maxSamples <- 50
    } else {
        maxSamples <- dots[["maxSamples"]]
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

    rowData <- rowData(se)
    colData <- colData(se)

    # Subset the tximport list
    txi <- bcbio(x, "tximport")
    txi[["abundance"]] <- txi[["abundance"]][genes, samples]
    txi[["counts"]] <- txi[["counts"]][genes, samples]
    txi[["length"]] <- txi[["length"]][genes, samples]

    # Obtain the counts from the updated tximport list
    raw <- txi[["counts"]]
    tmm <- .tmm(raw)
    tpm <- txi[["abundance"]]

    # Skip normalization option, for large datasets
    if (isTRUE(skipNorm)) {
        message("Skip re-normalization, just selecting samples and genes")
        # Only way to avoid disk space issue.
        # Direct subset of dds creates a huge file.
        dds <- .createDDS(txi, colData) %>%
            estimateSizeFactors()
        vst <- .subsetCounts(counts(x, "vst")[i, j], colData)
        rlog <- .subsetCounts(counts(x, "rlog")[i, j], colData)
        normalized <- counts(x, "normalized")[i, j]
    } else {
        # Fix for unexpected disk space issue (see constructor above)
        dds <- .createDDS(txi, colData)
        # DESeq2 will warn about empty design formula
        dds <- suppressWarnings(DESeq(dds))
        normalized <- counts(dds, normalized = TRUE)
    }

    # Update rlog and vst data
    if (nrow(colData) > maxSamples & !isTRUE(skipNorm)) {
        message("Many samples detected...skipping count transformations")
        rlog <- .subsetCounts(log2(tmm + 1), colData)
        vst <- .subsetCounts(log2(tmm + 1), colData)
    } else if (!isTRUE(skipNorm)) {
        message("Performing rlog transformation")
        rlog <- rlog(dds)
        message("Performing variance stabilizing transformation")
        vst <- varianceStabilizingTransformation(dds)
    }

    # Update featureCounts
    if (is.matrix(bcbio(x, "featureCounts"))) {
        featureCounts <- bcbio(x, "featureCounts") %>%
            .[genes, samples, drop = FALSE]
    } else {
        featureCounts <- NULL
    }

    # Update metadata ====
    metadata <- metadata(x)
    metadata[["metrics"]] <- metadata[["metrics"]] %>%
        .[.[["sampleID"]] %in% samples, , drop = FALSE]

    # Return `bcbioRNASeq`
    assays <- SimpleList(
        raw = raw,
        normalized = normalized,
        tpm = tpm,
        tmm = tmm,
        rlog = rlog,
        vst = vst)
    bcbio <- SimpleList(
        tximport = txi,
        DESeqDataSet = dds,
        featureCounts = featureCounts)
    new("bcbioRNASeq",
        SummarizedExperiment(
            assays = assays,
            colData = colData,
            rowData = rowData,
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
