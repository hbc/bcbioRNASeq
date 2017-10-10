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
#' plot53Bias(bcb, interestingGroup = "group")
#'
#' \dontrun{
#' # data.frame
#' metrics(bcb) %>%
#'     plot53Bias()
#' }
NULL



# Constructors ====
.plot53Bias <- function(
    object,
    interestingGroup = "sampleName",
    warnLimit = 2,
    flip = TRUE) {
    p <- ggplot(
        object,
        mapping = aes_string(
            x = "sampleName",
            y = "x53Bias",
            fill = interestingGroup)
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
setMethod("plot53Bias", "bcbioRNASeqANY", function(
    object,
    interestingGroup,
    warnLimit = 2,
    flip = TRUE) {
    if (is.null(metrics(object))) {
        return(NULL)
    }
    if (missing(interestingGroup)) {
        interestingGroup <- interestingGroups(object)[[1]]
    }
    .plot53Bias(
        metrics(object),
        interestingGroup = interestingGroup,
        warnLimit = warnLimit,
        flip = flip)
})



#' @rdname plot53Bias
#' @export
setMethod("plot53Bias", "data.frame", .plot53Bias)
