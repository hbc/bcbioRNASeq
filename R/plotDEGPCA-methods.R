# FIXME Define and recommend DESeqAnalysis method
# FIXME Include the alpha on the plot
# FIXME Check `DataFrame` return



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
#' plotDEGPCA(
#'     object = deseq_small@results,
#'     counts = deseq_small@transform,
#'     label = TRUE
#' )
NULL



#' @rdname plotDEGPCA
#' @export
setMethod(
    "plotDEGPCA",
    signature(
        object = "DESeqResults",
        counts = "DESeqTransform"
    ),
    function(
        object,
        counts,
        interestingGroups = NULL,
        alpha = NULL,
        lfcThreshold = 0L,
        direction = c("both", "up", "down"),
        color = getOption("bcbio.discrete.color", NULL),
        label = getOption("bcbio.label", FALSE),
        return = c("ggplot", "DataFrame")
    ) {
        validObject(object)
        validObject(counts)
        assert_are_identical(
            x = rownames(object),
            y = rownames(counts)
        )
        interestingGroups <- matchInterestingGroups(
            object = counts,
            interestingGroups = interestingGroups
        )
        interestingGroups(counts) <- interestingGroups
        if (is.null(alpha)) {
            alpha <- metadata(object)[["alpha"]]
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
            object,
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
            title = contrastName(object),
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
        object = "DESeqResults",
        counts = "bcbioRNASeq"
    ),
    function(
        object,
        counts,
        normalized = c("vst", "rlog", "tmm", "tpm", "rle"),
        ...
    ) {
        validObject(object)
        validObject(counts)
        assert_are_identical(
            x = rownames(object),
            y = rownames(counts)
        )
        normalized <- match.arg(normalized)
        message(paste("Using", normalized, "counts"))
        rse <- as(counts, "RangedSummarizedExperiment")
        assay(rse) <- counts(counts, normalized = normalized)
        plotDEGPCA(
            object = object,
            counts = rse,
            ...
        )
    }
)
