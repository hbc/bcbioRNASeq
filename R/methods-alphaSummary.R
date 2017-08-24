#' Print Summary Statistics of Alpha Level Cutoffs
#'
#' @rdname alphaSummary
#' @name alphaSummary
#'
#' @param alpha Numeric vector of desired alpha cutoffs.
#' @param contrast Character vector to use with [results()] function.
#' @param caption Character vector to add as caption to the table.
#' @return [kable].
#'
#' @examples
#' data(dds)
#' alphaSummary(dds)
NULL

.guessResults <- function(object, what, alpha){
    if (length(what) == 1)
        res <- results(object, name = what, alpha = alpha)
    else res <- results(object, contrast = what, alpha = alpha)
    res
}

# Methods ====
#' @rdname alphaSummary
#' @export
setMethod("alphaSummary", "DESeqDataSet", function(
    object,
    alpha = c(0.1, 0.05, 0.01, 1e-3, 1e-6),
    contrast = NULL,
    caption = NULL) {
    if (is.null(contrast))
        contrast <- resultsNames(object)[[2L]]
    if (is.null(caption))
        caption <- contrast
    lapply(seq_along(alpha), function(a) {
        .info <- capture.output(
            .guessResults(object, contrast, alpha[a]) %>%
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
        kable(., caption = paste(caption)) %>%
        show
})
