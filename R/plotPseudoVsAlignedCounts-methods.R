#' @name plotPseudoVsAlignedCounts
#' @author Michael Steinbaugh
#' @inherit AcidGenerics::plotPseudoVsAlignedCounts
#' @note Updated 2021-05-19.
#'
#' @note Currently supported for salmon or kallisto. The function will
#'   intentionally error for datasets containing aligned counts in the primary
#'   `counts` assay.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Passthrough to `AcidPlots::plotCountsCorrelationHeatmap()` when
#'   `genes = NULL` or `AcidPlots::plotCountsCorrelation()` when `genes` are
#'   defined.
#'
#' @examples
#' ## bcbioRNASeq ====
#' data(bcb)
#'
#' ## Correlation heatmap.
#' plotPseudoVsAlignedCounts(bcb)
#'
#' ## Individual genes.
#' ## Checking the most expressed aligned genes here.
#' aligned <- assay(bcb, i = "aligned")
#' genes <- names(tail(sort(rowSums(aligned)), n = 2L))
#' plotPseudoVsAlignedCounts(bcb, genes = genes)
NULL



## Updated 2020-01-20.
`plotPseudoVsAlignedCounts,bcbioRNASeq` <-  # nolint
    function(
        object,
        genes = NULL,
        title = "Pseudoaligned vs. aligned counts",
        ...
    ) {
        validObject(object)
        assert(
            .isGeneLevel(object),
            isSubset(c("counts", "aligned"), assayNames(object)),
            isString(title, nullOK = TRUE)
        )
        ## Coercing to SummarizedExperiment, for fast subsetting.
        object <- as(object, "RangedSummarizedExperiment")
        object <- humanize(object)
        ## Note that `i` here denotes the rows to keep.
        if (is.character(genes)) {
            assert(length(genes) <= 10L)
            i <- mapGenesToRownames(object = object, genes = genes)
        } else {
            ## Censor genes that aren't present in both, otherwise the
            ## correlation matrix calculation will fail.
            i <- !apply(
                X = assay(object, i = "aligned"),
                MARGIN = 1L,
                FUN = anyNA
            )
            if (sum(i) < nrow(object)) {
                n <- sum(!i, na.rm = TRUE)
                alertWarning(sprintf(
                    "Censoring %d %s containing an NA value.",
                    n,
                    ngettext(
                        n = n,
                        msg1 = "gene",
                        msg2 = "genes"
                    )
                ))
            }
        }
        object <- object[i, , drop = FALSE]
        pseudo <- assay(object, i = "counts")
        aligned <- assay(object, i = "aligned")
        assert(
            is.matrix(pseudo),
            !is.integer(pseudo),
            !anyNA(pseudo),
            is.matrix(aligned),
            is.integer(aligned),
            !anyNA(aligned)
        )
        if (is.character(genes)) {
            plotCountsCorrelation(
                x = pseudo,
                y = aligned,
                labels = list(
                    title = title
                ),
                .xname = "pseudoaligned",
                .yname = "aligned",
                ...
            )
        } else {
            plotCountsCorrelationHeatmap(
                x = pseudo,
                y = aligned,
                title = title,
                ...
            )
        }
    }



#' @rdname plotPseudoVsAlignedCounts
#' @export
setMethod(
    f = "plotPseudoVsAlignedCounts",
    signature = signature(object = "bcbioRNASeq"),
    definition = `plotPseudoVsAlignedCounts,bcbioRNASeq`
)
