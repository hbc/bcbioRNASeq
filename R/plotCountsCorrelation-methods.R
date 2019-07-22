#' Plot correlation of two count matrices
#'
#' @name plotCountsCorrelation
#' @export
#'
#' @inheritParams base::Extract
#' @inheritParams params
#'
#' @param ... Passthrough arguments to
#'   [`plotCorrelationHeatmap()`][acidplots::plotCorrelationHeatmap].
#'
#' @return `pheatmap` or `ggplot`.
#'   Returns `ggplot` when `genes` argument is defined, otherwise `pheatmap`.
#'
#' @examples
#' ## FIXME
NULL



# FIXME Allow user to subset using i, j.
`plotCountsCorrelation.matrix,matrix` <-  # nolint
    function(
        x,
        y,
        i = NULL,
        j = NULL,
        title = NULL,
        ...
    ) {
        validObject(x)
        validObject(y)
        assert(
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



#' @rdname plotCountsCorrelation
#' @export
setMethod(
    f = "plotCountsCorrelation",
    signature = signature(
        x = "matrix",
        y =  "matrix"
    ),
    definition = `plotCountsCorrelation.matrix,matrix`
)
