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
#' # DESeqResults, SummarizedExperiment ====
#' plotDEGPCA(
#'     results = res_small,
#'     counts = rld_small,
#'     label = TRUE
#' )
#'
#' # DESeqResults, bcbioRNASeq ====
#' plotDEGPCA(
#'     results = res_small,
#'     counts = bcb_small,
#'     label = TRUE
#' )
NULL



# Methods ======================================================================
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
        interestingGroups,
        lfc = 0L,
        color = scale_color_viridis(discrete = TRUE),
        label = FALSE,
        title = "deg pca",
        return = c("ggplot", "data.frame")
    ) {
        validObject(results)
        validObject(counts)
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(counts)
        }
        assertFormalInterestingGroups(colData(counts), interestingGroups)
        assert_is_a_number(lfc)
        assert_all_are_non_negative(lfc)
        assertIsColorScaleDiscreteOrNULL(color)
        assert_is_a_bool(label)
        assertIsAStringOrNULL(title)
        return <- match.arg(return)

        # Get the DE gene vector using `resultsTables()`
        list <- resultsTables(
            object = results,
            lfc = lfc,
            rowData = NULL,
            summary = FALSE,
            write = FALSE
        )
        deg <- c(
            list[["degLFCUp"]][["geneID"]],
            list[["degLFCDown"]][["geneID"]]
        )
        if (!length(deg)) {
            return(NULL)
        }

        # Subset the counts
        counts <- counts[deg, , drop = FALSE]

        # SummarizedExperiment method
        rse <- as(counts, "RangedSummarizedExperiment")
        plotPCA(
            object = rse,
            genes = rownames(rse),
            interestingGroups = interestingGroups,
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
        counts = "DESeqDataSet"
    ),
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
        normalized = c("rlog", "vst", "tmm", "tpm"),
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
)
