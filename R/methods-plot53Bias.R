#' Plot 5'->3' Bias
#'
#' @rdname plot53Bias
#' @author Michael Steinbaugh
#' @family Quality Control Plots
#' @inherit qcPlots
#'
#' @examples
#' data(bcb)
#'
#' # bcbioRNADataSet
#' plot53Bias(bcb)
#'
#' # data.frame
#' metrics <- metrics(bcb)
#' plot53Bias(metrics)



#' @rdname plot53Bias
.plot53Bias <- function(
    metrics,
    interestingGroup = "sampleName",
    warnLimit = 2L,
    flip = TRUE) {
    if (is.null(metrics)) return(NULL)
    p <- ggplot(metrics,
                aes_(x = ~sampleName,
                     y = ~x53Bias,
                     fill = as.name(interestingGroup))) +
        geom_bar(stat = "identity") +
        labs(title = "5'->3' bias",
             x = "sample",
             y = "5'->3' bias")
    if (!is.null(warnLimit)) {
        p <- p +
            geom_hline(alpha = qcLineAlpha,
                       color = qcWarnColor,
                       size = qcLineSize,
                       yintercept = warnLimit)
    }
    if (isTRUE(flip)) {
        p <- p + coord_flip()
    }
    p
}



#' @rdname plot53Bias
#' @export
setMethod("plot53Bias", "bcbioRNADataSet", .plot53Bias)



#' @rdname plot53Bias
#' @export
setMethod("plot53Bias", "data.frame", .plot53Bias)
