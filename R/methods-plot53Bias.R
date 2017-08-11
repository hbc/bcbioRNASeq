#' Plot 5'->3' Bias
#'
#' @rdname plot53Bias
#' @name plot53Bias
#'
#' @examples
#' data(bcb)
#'
#' # bcbioRNADataSet
#' plot53Bias(bcb)
#'
#' # data.frame
#' metrics(bcb) %>% plot53Bias
NULL



# Constructors ====
.plot53Bias <- function(
    object,
    interestingGroup = "sampleName",
    warnLimit = 2L,
    flip = TRUE) {
    if (is.null(object)) return(NULL)
    p <- ggplot(object,
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



# Methods ====
#' @rdname plot53Bias
#' @export
setMethod("plot53Bias", "bcbioRNADataSet", function(object, ...) {
    .plot53Bias(
        metrics(object),
        interestingGroup = .interestingGroup(object),
        ...)
})



#' @rdname plot53Bias
#' @export
setMethod("plot53Bias", "data.frame", .plot53Bias)
