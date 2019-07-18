#' Compare pseudoaligned counts to aligned counts.
#'
#' @note Currently supported for salmon or kallisto. The function will
#'   intentionally error for datasets containing aligned counts in the primary
#'   `counts` assay.
#'
#' @name plotPseudoVsAlignedCounts
#' @export
#'
#' @inheritParams params
#' @param genes `character`.
#'   Specific genes to plot. Currently limited to a max of 10 genes.
#' @param return `character(1)`.
#'   Return type. Uses [`match.arg()`][base::match.arg] internally. Only applies
#'   when `genes` argument is not defined.
#' @param ... Passthrough arguments.
#'   When `genes` argument is undefined, passes through to
#'   [`plotCorrelationHeatmap()`][acidplots::plotCorrelationHeatmap].
#'
#' @return `pheatmap` or `ggplot`.
#'   Returns `ggplot` when `genes` argument is defined, otherwise `pheatmap`.
#'
#' @examples
#' data(bcb)
#'
#' ## Correlation heatmap.
#' plotPseudoVsAlignedCounts(bcb)
#'
#' ## Individual genes.
#' ## Checking the most expressed aligned genes here.
#' genes <- assays(bcb) %>%
#'     .[["aligned"]] %>%
#'     rowSums() %>%
#'     sort() %>%
#'     tail(n = 2) %>%
#'     names()
#' plotPseudoVsAlignedCounts(bcb, genes = genes)
NULL



plotPseudoVsAlignedCounts.bcbioRNASeq <-  # nolint
    function(
        object,
        genes = NULL,
        title = "pseudo vs. aligned counts",
        return = c("plotCorrelationHeatmap", "cor"),
        ...
    ) {
        validObject(object)
        assert(
            isSubset(c("counts", "aligned"), assayNames(object)),
            isString(title, nullOK = TRUE)
        )
        return <- match.arg(return)

        pseudo <- assays(object)[["counts"]]
        aligned <- assays(object)[["aligned"]]

        assert(
            is.matrix(pseudo),
            !is.integer(pseudo),
            !anyNA(pseudo),
            is.matrix(aligned),
            is.integer(aligned)
        )

        # Censor genes that aren't present in both.
        censor <- apply(X = aligned, MARGIN = 1L, FUN = anyNA)
        message(sprintf("Censoring %d genes.", sum(censor)))
        keep <- !censor
        pseudo <- pseudo[keep, , drop = FALSE]
        aligned <- aligned[keep, , drop = FALSE]
        assert(!anyNA(aligned))

        if (is.character(genes)) {
            assert(
                isCharacter(genes),
                length(genes) <= 10L,
                isSubset(genes, rownames(pseudo)),
                isSubset(genes, rownames(aligned))
            )
            pseudo %<>%
                .[genes, , drop = FALSE] %>%
                meltCounts(minCounts = NULL) %>%
                ungroup() %>%
                mutate_if(is.factor, as.character) %>%
                mutate(!!sym("type") := "pseudo")
            aligned %<>%
                .[genes, , drop = FALSE] %>%
                meltCounts(minCounts = NULL) %>%
                ungroup() %>%
                mutate_if(is.factor, as.character) %>%
                mutate(!!sym("type") := "aligned")
            data <- bind_rows(pseudo, aligned)
            ggplot(
                data = data,
                mapping = aes(
                    x = !!sym("colname"),
                    y = !!sym("counts"),
                    colour = !!sym("type")
                )
            ) +
                geom_point() +
                facet_wrap(facets = sym("rowname"), scales = "free_y") +
                labs(x = "sample", title = title)
        } else {
            # Correlation heatmap.
            # Consider limiting to the top n genes expressed instead.
            cor <- cor(x = pseudo, y = aligned, method = "pearson")
            assert(!anyNA(cor))
            se <- SummarizedExperiment(
                assays = list(cor = cor),
                colData = colData(object)
            )
            plotCorrelationHeatmap(
                object = se,
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
    definition = plotPseudoVsAlignedCounts.bcbioRNASeq
)
