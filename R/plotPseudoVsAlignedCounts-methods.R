## FIXME Switch to labels.



#' Compare pseudoaligned counts to aligned counts.
#'
#' @name plotPseudoVsAlignedCounts
#' @note Currently supported for salmon or kallisto. The function will
#'   intentionally error for datasets containing aligned counts in the primary
#'   `counts` assay.
#' @note Updated 2019-09-15.
#'
#' @inheritParams acidroxygen::params
#' @param ... Passthrough to [acidplots::plotCountsCorrelationHeatmap()] when
#'   `genes = NULL` or [acidplots::plotCountsCorrelation()] when `genes` are
#'   defined.
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



## Updated 2019-09-15.
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
        pseudo <- assays(object)[["counts"]]
        aligned <- assays(object)[["aligned"]]
        assert(
            is.matrix(pseudo),
            !is.integer(pseudo),
            !anyNA(pseudo),
            is.matrix(aligned),
            is.integer(aligned)
        )
        if (is.character(genes)) {
            assert(
                isCharacter(genes),
                length(genes) <= 10L,
                isSubset(genes, rownames(pseudo)),
                isSubset(genes, rownames(aligned))
            )
            pseudo <- pseudo[genes, , drop = FALSE]
            aligned <- aligned[genes, , drop = FALSE]
            ## FIXME Title isn't supported here yet.
            plotCountsCorrelation(
                x = pseudo,
                y = aligned,
                title = title,
                ...
            )
        } else {
            ## Censor genes that aren't present in both, otherwise the
            ## correlation matrix calculation will fail.
            censor <- apply(X = aligned, MARGIN = 1L, FUN = anyNA)
            message(sprintf("Censoring %d genes.", sum(censor)))
            keep <- !censor
            pseudo <- pseudo[keep, , drop = FALSE]
            aligned <- aligned[keep, , drop = FALSE]
            assert(!anyNA(aligned))
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
