#' Plot 5'->3' Bias
#'
#' @rdname plot53Bias
#' @name plot53Bias
#' @family Quality Control Plots
#' @author Michael Steinbaugh
#'
#' @inherit qcPlots
#'
#' @inheritParams AllGenerics
#'
#' @examples
#' data(bcb)
#'
#' # bcbioRNASeq
#' plot53Bias(bcb)
#'
#' \dontrun{
#' plot53Bias(bcb, interestingGroups = "group")
#'
#' # data.frame
#' metrics(bcb) %>%
#'     plot53Bias()
#' }
NULL



# Constructors ====
#' @importFrom viridis scale_fill_viridis
.plot53Bias <- function(
    object,
    interestingGroups = "sampleName",
    warnLimit = 2,
    flip = TRUE) {
    p <- ggplot(
        object,
        mapping = aes_string(
            x = "sampleName",
            y = "x53Bias",
            fill = interestingGroups)
    ) +
        geom_bar(stat = "identity") +
        labs(title = "5'->3' bias",
             x = "sample",
             y = "5'->3' bias") +
        scale_fill_viridis(discrete = TRUE)
    if (!is.null(warnLimit)) {
        p <- p +
            qcWarnLine(warnLimit)
    }
    if (isTRUE(flip)) {
        p <- p +
            coord_flip()
    }
    p
}



# Methods ====
#' @rdname plot53Bias
#' @export
setMethod(
    "plot53Bias",
    signature("bcbioRNASeqANY"),
    function(
        object,
        interestingGroups,
        warnLimit = 2,
        flip = TRUE) {
        if (is.null(metrics(object))) {
            return(NULL)
        }
        if (missing(interestingGroups)) {
            interestingGroups <-
                metadata(object)[["interestingGroups"]][[1]]
        }
        .plot53Bias(
            metrics(object),
            interestingGroups = interestingGroups,
            warnLimit = warnLimit,
            flip = flip)
    })



#' @rdname plot53Bias
#' @export
setMethod(
    "plot53Bias",
    signature("data.frame"),
    .plot53Bias)
