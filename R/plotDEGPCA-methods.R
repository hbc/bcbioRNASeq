# FIXME Include the alpha on the plot



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
NULL



#' @rdname plotDEGPCA
#' @export
setMethod(
    "plotDEGPCA",
    signature(
        results = "DESeqResults",
        counts = "SummarizedExperiment"
    ),
    function(
        results,
        counts,
        interestingGroups = NULL,
        alpha = NULL,
        lfcThreshold = 0L,
        direction = c("both", "up", "down"),
        color = getOption("bcbio.discrete.color", NULL),
        label = getOption("bcbio.label", FALSE),
        return = c("ggplot", "data.frame")
    ) {
        validObject(results)
        validObject(counts)
        assert_are_identical(
            x = rownames(results),
            y = rownames(counts)
        )
        interestingGroups <- matchInterestingGroups(
            object = counts,
            interestingGroups = interestingGroups
        )
        interestingGroups(counts) <- interestingGroups
        if (is.null(alpha)) {
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
            warning("No significant DEGs to plot", call. = FALSE)
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
)



#' @rdname plotDEGPCA
#' @export
setMethod(
    "plotDEGPCA",
    signature(
        results = "DESeqResults",
        counts = "bcbioRNASeq"
    ),
    function(
        results,
        counts,
        normalized = c("vst", "rlog", "tmm", "tpm", "rle"),
        ...
    ) {
        validObject(results)
        validObject(counts)
        assert_are_identical(
            x = rownames(results),
            y = rownames(counts)
        )
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
)



#' @rdname plotDEGPCA
#' @export
setMethod(
    "plotDEGPCA",
    signature(
        results = "DESeqResults",
        counts = "DESeqDataSet"
    ),
    function(
        results,
        counts,
        ...
    ) {
        validObject(results)
        validObject(counts)
        assert_are_identical(
            x = rownames(results),
            y = rownames(counts)
        )
        message("Using normalized counts")
        rse <- as(counts, "RangedSummarizedExperiment")
        assay(rse) <- counts(counts, normalized = TRUE)
        plotDEGPCA(
            results = results,
            counts = rse,
            ...
        )
    }
)



#' @rdname plotDEGPCA
#' @export
setMethod(
    "plotDEGPCA",
    signature(
        results = "DESeqResults",
        counts = "DESeqTransform"
    ),
    getMethod(
        "plotDEGPCA",
        signature(
            results = "DESeqResults",
            counts = "SummarizedExperiment"
        )
    )
)
