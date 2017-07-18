#' Make groups of genes using expression profile
#'
#' [DEGreport::degPatterns()] wrapper supporting a [bcbioRNADataSet].
#'
#' @rdname plot_pattern
#' @author Lorena Pantano
#'
#' @param bcb [bcbioRNADataSet] object.
#' @param res Table with padj as column and gene names as [rownames].
#' @param fdr Float cutoff to consider genes significant.
#' @param n Integer maximum number of genes to consider.
#' @param ... Additional parameters, passed to [DEGreport::degPatterns()].
#'
#' @return [ggplot].
#' @export
setMethod("plot_pattern", "bcbioRNADataSet", function(
    object,
    res,
    fdr = 0.1,
    n = NULL,
    ...) {
    res <- as.data.frame(res)
    .order <- res[order(res[["padj"]]), ]
    sign <- row.names(.order)[!is.na(.order[["padj"]])]
    if (!is.null(n)) {
        sign <- sign[1L:n]
    }
    degPatterns(counts(object, normalized = "rlog")[sign, ],
                metadata = colData(object),
                ...)
})
