#' Plot MA
#'
#' Wrapper function for \code{DESeq2::plotMA()} that enforces a symmetric y-axis
#' and an automatic title
#'
#' @author Michael Steinbaugh
#'
#' @import DESeq2
#'
#' @param res \code{DESeqResults} object
#' @param ylim Y-axis maximum (single integer)
#'
#' @return MA plot
#' @export
plot_ma <- function(res, ylim = 2) {
    if (class(res)[1] != "DESeqResults") {
        stop("DESeqResults required")
    }
    DESeq2::plotMA(res,
                   main = res_contrast_name(res),
                   ylim = c(-ylim, ylim))
}
