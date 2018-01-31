#' Bracket-Based Subsetting
#'
#' Extract genes by row and samples by column from a [bcbioRNASeq] object. The
#' internal [DESeqDataSet] and count transformations are rescaled automatically.
#' DESeq2 transformations can be disabled on large subset operations by setting
#' `transform = FALSE`.
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
#'     file.path("extdata", "bcb.rda"),
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
#' # Skip DESeq2 transformations
#' bcb[ensgene, samples, transform = FALSE]
NULL



# Constructors =================================================================
#' @importFrom DESeq2 DESeq estimateSizeFactors rlog
#'   varianceStabilizingTransformation
#' @importFrom dplyr mutate_if
#' @importFrom S4Vectors metadata SimpleList
#' @importFrom tibble column_to_rownames rownames_to_column
.subset <- function(x, i, j, ..., drop = FALSE) {
    # Genes (rows)
    if (missing(i)) {
        i <- 1L:nrow(x)
    }
    if (length(i) < 2L) {
        abort("At least 2 genes are required")
    }

    # Samples (columns)
    if (missing(j)) {
        j <- 1L:ncol(x)
    }
    if (length(j) < 2L) {
        abort("At least 2 samples are required")
    }

    # Early return if dimensions are unmodified
    if (identical(dim(x), c(length(i), length(j)))) return(x)

    dots <- list(...)
    if (!identical(dots[["transform"]], FALSE)) {
        transform <- TRUE
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
    inform("Updating internal DESeqDataSet")
    dds <- bcbio(x, "DESeqDataSet")
    dds <- dds[genes, samples]
    colData(dds) <- colData
    # Skip normalization option, for large datasets
    if (!isTRUE(transform)) {
        inform(paste(
            "Skipping DESeq2 transformations",
            "just selecting samples and genes"
        ))
        dds <- estimateSizeFactors(dds)
        vst <- NULL
        rlog <- NULL
    } else {
        # DESeq2 will warn about empty design formula, if set
        dds <- suppressWarnings(DESeq(dds))
        inform("Performing rlog transformation")
        rlog <- rlog(dds)
        inform("Performing variance stabilizing transformation")
        vst <- varianceStabilizingTransformation(dds)
    }

    # featureCounts
    featureCounts <- bcbio(x, "featureCounts")
    if (is.matrix(featureCounts)) {
        # Genes detected by pseudo-alignment (e.g. salmon) will differ from
        # the genes detected by alignment with STAR/featureCounts. Need to
        # obtain the intersect here to avoid a subset error. Note that these
        # counts are only used to generate the quality control metrics for
        # MultiQC and should not be used for analysis.
        commonGenes <- intersect(rownames(featureCounts), genes)
        featureCounts <- featureCounts[commonGenes, samples, drop = FALSE]
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



# Methods ======================================================================
#' @rdname subset
#' @export
setMethod(
    "[",
    signature(
        x = "bcbioRNASeq",
        i = "ANY",
        j = "ANY"),
    function(x, i, j, ..., drop = FALSE) {
        .subset(x, i, j, ..., drop)
    })
