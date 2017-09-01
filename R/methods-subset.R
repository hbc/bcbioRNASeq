#' Bracket-Based Subsetting
#'
#' Extract genes by row and samples by column from a [bcbioRNADataSet]. The
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
#' @return [bcbioRNADataSet].
#'
#' @seealso `help("[", "base")`.
#'
#' @examples
#' data(bcb)
#' genes <- 1L:50L
#' samples <- c("group1_1", "group1_2")
#'
#' # Subset by sample name
#' bcb[, samples]
#'
#' # Subset by gene list
#' bcb[genes, ]
#'
#' # Subset by both genes and samples
#' bcb[genes, samples]
NULL



# Constructors ====
# This operation must be placed outside of the S4 method dispatch. Otherwise,
# the resulting subset object will be ~2X the expected size on disk when saving,
# for an unknown reason.
.createDDS <- function(txi, tmpData) {
    DESeqDataSetFromTximport(
        txi = txi,
        colData = tmpData,
        design = formula(~1L))
}

.countSubset <- function(x, tmpData){
        DESeqTransform(
            SummarizedExperiment(assays = x,
                                 colData = tmpData))
}

# Methods ====
#' @rdname subset
#' @export
setMethod(
    "[",
    signature(x = "bcbioRNADataSet", i = "ANY", j = "ANY"),
    function(x, i, j, ..., drop = FALSE) {
        if (missing(i)) {
            i <- 1L:nrow(x)
        }
        if (missing(j)) {
            j <- 1L:ncol(x)
        }
        dots <- list(...)
        if (is.null(dots[["maxSamples"]])) {
            maxSamples <- 50L
        } else {
            maxSamples <- dots[["maxSamples"]]
        }

        if (is.null(dots[["skipNorm"]])) {
            skipNorm <- FALSE
        } else {
            skipNorm <- dots[["skipNorm"]]
        }

        # Subset SE object ====
        se <- SummarizedExperiment(
            assays = SimpleList(counts(x)),
            rowData = rowData(x),
            colData = colData(x),
            metadata = metadata(x))

        tmp <- se[i, j, drop = drop]
        tmpGenes <- row.names(tmp)
        tmpRow <- rowData(tmp) %>% as.data.frame
        tmpData <- colData(tmp) %>% as.data.frame
        tmpTxi <- bcbio(x, "tximport")

        # Subset tximport data ====
        txi <- SimpleList(
            abundance = tmpTxi[["abundance"]] %>%
                .[tmpGenes, tmpData[["sampleID"]]],
            counts = tmpTxi[["counts"]] %>%
                .[tmpGenes, tmpData[["sampleID"]]],
            length = tmpTxi[["length"]] %>%
                .[tmpGenes, tmpData[["sampleID"]]],
            countsFromAbundance = tmpTxi[["countsFromAbundance"]]
        )
        rawCounts <- txi[["counts"]]
        tmm <- .tmm(rawCounts)
        tpm <- txi[["abundance"]]

        if (skipNorm){
            message("Skip re-normalization, just selecting samples and genes.")
            # To Fix if we find the solution.
            # Only way to avoid disk space issue.
            # Direct subset of dds create a huge file.
            dds <- .createDDS(txi, tmpData)  %>%
                estimateSizeFactors()
            vst <- .countSubset(counts(x, "vst")[i, j], tmpData)
            rlog <- .countSubset(counts(x, "rlog")[i, j], tmpData)
            normalizedCounts <- counts(x, "normalized")[i, j]
        } else {
            # Fix for unexpected disk space issue (see constructor above)
            dds <- .createDDS(txi, tmpData)  %>%
                DESeq
            normalizedCounts <- counts(dds, normalized = TRUE)
        }

        # rlog & variance ====
        if (nrow(tmpData) > maxSamples & !skipNorm) {
            message("Many samples detected...skipping count transformations")
            rlog <- .countSubset(log2(tmm + 1L), tmpData)
            vst <- .countSubset(log2(tmm + 1L), tmpData)
        } else if (!skipNorm) {
            message("Performing rlog transformation")
            rlog <- rlog(dds)
            message("Performing variance stabilizing transformation")
            vst <- varianceStabilizingTransformation(dds)
        }

        if (is.matrix(bcbio(x, "featureCounts"))) {
            tmpFC <- bcbio(x, "featureCounts") %>%
                .[tmpGenes, tmpData[["sampleID"]]]
        } else {
            tmpFC <- bcbio(x, "featureCounts")
        }

        # Subset Metrics ====
        tmpMetadata <- metadata(x)
        tmpMetrics <- tmpMetadata[["metrics"]] %>%
            .[.[["sampleID"]] %in% tmpData[["sampleID"]], ]
        tmpMetadata[["metrics"]] <- tmpMetrics

        tmpAssays <- SimpleList(
            raw = rawCounts,
            normalized = normalizedCounts,
            tpm = tpm,
            tmm = tmm,
            rlog = rlog,
            vst = vst)

        # Slot additional callers
        extra <- list(tximport = txi,
                      DESeqDataSet = dds,
                      featureCounts = tmpFC)

        # bcbioRNADataSet ====
        # FIXME Calling an internal function with `:::` isn't valid
        bcb <- BiocGenerics:::replaceSlots(
            x,
            assays = tmpAssays,
            colData = tmpData,
            elementMetadata = tmpRow,
            callers = extra,
            metadata = tmpMetadata,
            check = FALSE)
        bcb
    })
