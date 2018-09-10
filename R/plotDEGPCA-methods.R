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
#'     object = deseq_small@lfcShrink[[1L]],
#'     counts = deseq_small@transform,
#'     label = TRUE
#' )
NULL



# FIXME Use the formals defined in `plotPCA`...
.plotDEGPCA.DESeqResults.DESeqTransform <-  # nolint
    function(
        object,
        counts,
        interestingGroups = NULL,
        direction = c("both", "up", "down"),
        color = getOption("bcbio.discrete.color", NULL),
        label = getOption("bcbio.label", TRUE),
        return = c("ggplot", "DataFrame")
    ) {
        assert_is_all_of(object, "DESeqResults")
        validObject(object)
        assert_is_all_of(counts, "DESeqTransform")
        validObject(counts)
        assert_are_identical(rownames(object), rownames(counts))
        interestingGroups <- matchInterestingGroups(
            object = counts,
            interestingGroups = interestingGroups
        )
        interestingGroups(counts) <- interestingGroups
        alpha <- metadata(object)[["alpha"]]
        assert_is_a_number(alpha)
        lfcThreshold <- metadata(object)[["lfcThreshold"]]
        assert_is_a_number(lfcThreshold)
        assert_all_are_non_negative(lfcThreshold)
        assertIsColorScaleDiscreteOrNULL(color)
        direction <- match.arg(direction)
        assert_is_a_bool(label)
        return <- match.arg(return)

        # Get DEG vector using DEGreport.
        if (direction == "both") {
            direction <- NULL
        }
        deg <- significants(
            object,
            padj = alpha,
            fc = lfcThreshold,
            direction = direction
        )

        # Early return if there are no DEGs.
        if (!length(deg) > 0L) {
            warning("No significant DEGs to plot", call. = FALSE)
            return(invisible())
        }

        # Using our internal `SummarizedExperiment` method.
        do.call(
            what = .plotPCA.SummarizedExperiment,
            args = list(
                object = counts[deg, , drop = FALSE],
                interestingGroups = interestingGroups,
                ntop = Inf,
                label = label,
                title = contrastName(object),
                subtitle = paste(length(deg), "genes"),
                return = return
            )
        )
    }



.plotDEGPCA.DESeqResults.bcbioRNASeq <-  # nolint
    function(
        object,
        counts,
        normalized = c("vst", "rlog", "tmm", "tpm", "rle"),
        ...
    ) {
        validObject(object)
        validObject(counts)
        assert_are_identical(rownames(object), rownames(counts))
        normalized <- match.arg(normalized)
        message(paste("Using", normalized, "counts"))
        rse <- as(counts, "RangedSummarizedExperiment")
        assays(rse) <- list(counts(counts, normalized = normalized))
        # Handing off to `DESeqResults,DESeqTransform` method.
        do.call(
            what = plotDEGPCA,
            args = list(
                object = object,
                counts = DESeqTransform(rse),
                ...
            )
        )
    }
# FIXME Set the formals.



.plotDEGPCA.DESeqAnalysis <-  # nolint
    function(
        object,
        counts = NULL,
        results = 1L
    ) {
        results <- object@lfcShrink[[results]]
        counts <- object@transform
        args <- setArgsToDoCall(
            args = list(
                object = results,
                counts = counts
            ),
            removeArgs = "results",
            call = matchCall()
        )
        do.call(what = plotDEGPCA, args = args)
    }



#' @rdname plotDEGPCA
#' @export
setMethod(
    f = "plotDEGPCA",
    signature = signature(
        object = "DESeqAnalysis",
        counts = "missingOrNULL"
    ),
    definition = .plotDEGPCA.DESeqAnalysis
)



#' @rdname plotDEGPCA
#' @export
setMethod(
    f = "plotDEGPCA",
    signature = signature(
        object = "DESeqResults",
        counts = "DESeqTransform"
    ),
    definition = .plotDEGPCA.DESeqResults.DESeqTransform
)



#' @rdname plotDEGPCA
#' @export
setMethod(
    f = "plotDEGPCA",
    signature = signature(
        object = "DESeqResults",
        counts = "bcbioRNASeq"
    ),
    definition = .plotDEGPCA.DESeqResults.bcbioRNASeq
)
