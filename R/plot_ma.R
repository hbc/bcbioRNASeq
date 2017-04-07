#' Differential expression plots
#'
#' @rdname de_plots
#' @author Michael Steinbaugh
#'
#' @description Wrapper for \code{DESeq2::plotMA()} that enforces a symmetric
#'   y-axis and an automatic title
#'
#' @import DESeq2
#'
#' @param res \code{DESeqResults} object
#' @param ylim Y-axis maximum (single integer)
#'
#' @return MA plot
#' @export
plot_ma <- function(res, ylim = 2) {
    check_res(res)

    name <- deparse(substitute(res))
    contrast_name <- res_contrast_name(res)

    plot <- DESeq2::plotMA(
        res,
        main = paste0(name, ": ", contrast_name),
        ylim = c(-ylim, ylim))

    show(plot)
}
