#' Plot DEG PCA
#'
#' @name plotDEGPCA
#' @family Differential Expression Functions
#' @author Michael Steinbaugh
#'
#' @inherit plotPCA
#' @inheritParams general
#'
#' @examples
#' # DESeqResults, DESeqTransform ====
#' vst_small <- DESeq2::varianceStabilizingTransformation(dds_small)
#' plotDEGPCA(
#'     results = res_small,
#'     counts = vst_small,
#'     label = TRUE
#' )
#'
#' # DESeqResults, bcbioRNASeq ====
#' plotDEGPCA(
#'     results = res_small,
#'     counts = bcb_small,
#'     normalized = "vst",
#'     label = TRUE
#' )
NULL



#' @rdname plotDEGPCA
#' @export
setGeneric(
    "plotDEGPCA",
    function(results, counts, ...) {
        standardGeneric("plotDEGPCA")
    }
)



`plotDEGPCA.DESeqResults,SummarizedExperiment` <-  # nolint
    function(
        results,
        counts,
        interestingGroups,
        alpha,
        lfcThreshold = 0L,
        direction = c("both", "up", "down"),
        color = getOption("bcbio.discrete.color", NULL),
        label = getOption("bcbio.label", FALSE),
        return = c("ggplot", "data.frame")
    ) {
        validObject(results)
        validObject(counts)
        interestingGroups <- matchInterestingGroups(
            object = counts,
            interestingGroups = interestingGroups
        )
        if (missing(alpha)) {
            alpha <- metadata(results)[["alpha"]]
        }
        assert_is_a_number(alpha)
        assert_is_a_number(lfcThreshold)
        assert_all_are_non_negative(lfcThreshold)
        assertIsColorScaleDiscreteOrNULL(color)
        direction <- match.arg(direction)
        assert_is_a_bool(label)
        return <- match.arg(return)

        # Get DEG vector using DEGreport
        if (direction == "both") {
            direction <- NULL
        }
        deg <- significants(
            results,
            padj = alpha,
            fc = lfcThreshold,
            direction = direction
        )

        # Early return if there are no DEGs
        if (!length(deg)) {
            warning("No significant DEGs to plot")
            return(invisible())
        }

        # Subset the counts
        counts <- counts[deg, , drop = FALSE]

        # SummarizedExperiment method
        rse <- as(counts, "RangedSummarizedExperiment")
        plotPCA(
            object = rse,
            interestingGroups = interestingGroups,
            ntop = Inf,
            label = label,
            title = contrastName(results),
            subtitle = paste(nrow(rse), "genes"),
            return = return
        )
    }



#' @rdname plotDEGPCA
#' @export
setMethod(
    f = "plotDEGPCA",
    signature = signature(
        results = "DESeqResults",
        counts = "SummarizedExperiment"
    ),
    definition = `plotDEGPCA.DESeqResults,SummarizedExperiment`
)



`plotDEGPCA.DESeqResults,bcbioRNASeq` <-  # nolint
    function(
        results,
        counts,
        normalized = c("vst", "rlog", "tmm", "tpm", "rle"),
        ...
    ) {
        validObject(counts)
        normalized <- match.arg(normalized)
        message(paste("Using", normalized, "counts"))
        rse <- as(counts, "RangedSummarizedExperiment")
        assay(rse) <- counts(counts, normalized = normalized)
        plotDEGPCA(
            results = results,
            counts = rse,
            ...
        )
    }



#' @rdname plotDEGPCA
#' @export
setMethod(
    f = "plotDEGPCA",
    signature = signature(
        results = "DESeqResults",
        counts = "bcbioRNASeq"
    ),
    definition = `plotDEGPCA.DESeqResults,bcbioRNASeq`
)



`plotDEGPCA.DESeqResults,DESeqDataSet` <-  # nolint
    function(
        results,
        counts,
        ...
    ) {
        validObject(counts)
        message("Using normalized counts")
        rse <- as(counts, "RangedSummarizedExperiment")
        assay(rse) <- counts(counts, normalized = TRUE)
        plotDEGPCA(
            results = results,
            counts = rse,
            ...
        )
    }



#' @rdname plotDEGPCA
#' @export
setMethod(
    f = "plotDEGPCA",
    signature = signature(
        results = "DESeqResults",
        counts = "DESeqDataSet"
    ),
    definition = `plotDEGPCA.DESeqResults,DESeqDataSet`
)



`plotDEGPCA.DESeqResults,DESeqTransform` <-  # nolint
    `plotDEGPCA.DESeqResults,SummarizedExperiment`



#' @rdname plotDEGPCA
#' @export
setMethod(
    f = "plotDEGPCA",
    signature = signature(
        results = "DESeqResults",
        counts = "DESeqTransform"
    ),
    definition = `plotDEGPCA.DESeqResults,DESeqTransform`
)
