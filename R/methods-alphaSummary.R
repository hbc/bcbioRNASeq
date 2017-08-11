#' Print Summary Statistics of Alpha Level Cutoffs
#'
#' @rdname alphaSummary
#' @name alphaSummary
#'
#' @param alpha Numeric vector of desired alpha cutoffs.
#' @param contrast Character vector to use with [results()] function.
#'
#' @return [kable].
#'
#' @examples
#' data(dds)
#' alphaSummary(dds)
NULL



# Methods ====
#' @rdname alphaSummary
#' @export
setMethod("alphaSummary", "DESeqDataSet", function(
    object,
    alpha = c(0.1, 0.05, 0.01, 1e-3, 1e-6),
    contrast = NULL) {
    if (is.null(contrast))
        contrast <- strsplit(resultsNames(object)[[2L]], "_") %>%
            unlist %>%
            .[c(1L, 2L, 4L)]
    lapply(seq_along(alpha), function(a) {
        .info <- capture.output(
            results(object, contrast = contrast, alpha = alpha[a]) %>%
                summary)[4L:8L]
        .parse <- sapply(
            .info, function(i) {
                gsub(".*:", "", i)
            })[1L:4L]
        .parse <- c(.parse, .info[[5L]])
        data.frame(alpha = as.vector(.parse))
    }) %>%
        bind_cols %>%
        set_colnames(alpha) %>%
        set_rownames(c("LFC > 0 (up)",
                       "LFC < 0 (down)",
                       "outliers",
                       "low counts",
                       "cutoff")) %>%
        kable %>%
        show
})
