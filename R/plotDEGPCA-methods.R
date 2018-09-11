#' Plot DEG PCA
#'
#' @name plotDEGPCA
#' @family Differential Expression Functions
#' @author Michael Steinbaugh
#' @include plotPCA-methods.R
#'
#' @inherit plotPCA
#' @inheritParams general
#'
#' @examples
#' # DESeqAnalysis ====
#' plotDEGPCA(deseq_small)
#'
#' # DESeqResults, DESeqTransform ====
#' plotDEGPCA(
#'     object = deseq_small@lfcShrink[[1L]],
#'     counts = deseq_small@transform
#' )
#'
#' # DESeqResults, bcbioRNASeq ====
#' plotDEGPCA(
#'     object = deseq_small@lfcShrink[[1L]],
#'     counts = bcb_small,
#'     normalized = "vst"
#' )
NULL



.plotDEGPCA.DESeqResults.DESeqTransform <-  # nolint
    function(
        object,
        counts,
        direction = c("both", "up", "down")
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
                subtitle = paste(
                    length(deg), "genes;",
                    "alpha <", alpha
                ),
                return = return
            )
        )
    }

# Assign the formals.
f1 <- formals(.plotDEGPCA.DESeqResults.DESeqTransform)
f2 <- formals(.plotPCA.SummarizedExperiment)
f2 <- f2[c("interestingGroups", "color", "label", "return")]
f <- c(f1, f2)
formals(.plotDEGPCA.DESeqResults.DESeqTransform) <- f



.plotDEGPCA.DESeqResults.bcbioRNASeq <-  # nolint
    function(
        object,
        counts,
        normalized = c("vst", "rlog", "tmm", "tpm", "rle")
    ) {
        validObject(object)
        validObject(counts)
        assert_are_identical(rownames(object), rownames(counts))
        normalized <- match.arg(normalized)
        message(paste("Using", normalized, "counts"))
        rse <- as(counts, "RangedSummarizedExperiment")
        assays(rse) <- list(counts(counts, normalized = normalized))
        # Handing off to `DESeqResults,DESeqTransform` method.
        args <- setArgsToDoCall(
            args = list(
                object = object,
                counts = DESeqTransform(rse)
            ),
            removeArgs = "normalized",
            call = matchCall()
        )
        do.call(what = plotDEGPCA, args = args)
    }

# Assign the formals.
f1 <- formals(.plotDEGPCA.DESeqResults.bcbioRNASeq)
f2 <- formals(.plotDEGPCA.DESeqResults.DESeqTransform)
f2 <- f2[setdiff(names(f2), names(f1))]
f <- c(f1, f2)
formals(.plotDEGPCA.DESeqResults.bcbioRNASeq) <- f



.plotDEGPCA.DESeqAnalysis <-  # nolint
    function(
        object,
        counts = NULL,
        results
    ) {
        results <- .matchResults(
            object = object,
            results = results
        )
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

# Assign the formals.
f1 <- formals(.plotDEGPCA.DESeqAnalysis)
f2 <- formals(.plotDEGPCA.DESeqResults.DESeqTransform)
f2 <- f2[setdiff(names(f2), names(f1))]
f <- c(f1, f2)
formals(.plotDEGPCA.DESeqAnalysis) <- f


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
