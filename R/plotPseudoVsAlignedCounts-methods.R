#' Compare pseudoaligned counts to aligned counts.
#'
#' @name plotPseudoVsAlignedCounts
#' @note Currently supported for salmon or kallisto. The function will
#'   intentionally error for datasets containing aligned counts in the primary
#'   `counts` assay.
#' @note Updated 2020-01-20.
#'
#' @inheritParams acidroxygen::params
#' @param ... Passthrough to [acidplots::plotCountsCorrelationHeatmap()] when
#'   `genes = NULL` or [acidplots::plotCountsCorrelation()] when `genes` are
#'   defined.
#'
#' @return Plot.
#'
#' @examples
#' data(bcb)
#'
#' ## Correlation heatmap.
#' plotPseudoVsAlignedCounts(bcb)
#'
#' ## Individual genes.
#' ## Checking the most expressed aligned genes here.
#' aligned <- assay(bcb, i = "aligned")
#' genes <- rowSums(aligned) %>% sort %>% tail(n = 2L) %>% names()
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
                cli_alert_warning(sprintf(
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
    signature = signature("bcbioRNASeq"),
    definition = `plotPseudoVsAlignedCounts,bcbioRNASeq`
)
