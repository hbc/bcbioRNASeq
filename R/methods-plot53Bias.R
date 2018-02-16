#' Plot 5'->3' Bias
#'
#' @rdname plot53Bias
#' @name plot53Bias
#' @family Quality Control Plots
#' @author Michael Steinbaugh
#'
#' @inherit plotTotalReads
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq
#' plot53Bias(bcb)
#' plot53Bias(
#'     bcb,
#'     interestingGroups = "sampleName",
#'     fill = NULL)
#'
#' # data.frame
#' df <- metrics(bcb)
#' plot53Bias(df)
NULL



# Constructors =================================================================
#' @importFrom bcbioBase uniteInterestingGroups
#' @importFrom dplyr rename
#' @importFrom ggplot2 aes_string coord_flip geom_bar ggplot guides labs
.plot53Bias <- function(
    object,
    interestingGroups = "sampleName",
    warnLimit = 2L,
    fill = viridis::scale_fill_viridis(discrete = TRUE),
    flip = TRUE,
    title = TRUE) {
    assert_is_data.frame(object)
    assert_is_character(interestingGroups)
    assert_is_an_implicit_integer(warnLimit)
    assert_is_any_of(fill, c("ScaleDiscrete", NULL))
    assert_is_a_bool(flip)
    assert_is_any_of(title, c("character", "NULL"))

    data <- uniteInterestingGroups(object, interestingGroups)

    # Legacy code: make sure `x53Bias` is camel sanitized to `x5x3Bias`.
    # The internal camel method has been updated in basejump 0.1.11.
    if ("x53Bias" %in% colnames(data)) {
        data <- rename(data, "x5x3Bias" = "x53Bias")
    }

    ylab <- "5'->3' bias"
    if (isTRUE(title)) {
        title <- ylab
    } else if (!is_a_string(title)) {
        title <- NULL
    }

    p <- ggplot(
        data = data,
        mapping = aes_string(
            x = "sampleName",
            y = "x5x3Bias",
            fill = "interestingGroups")
    ) +
        geom_bar(stat = "identity") +
        labs(
            title = title,
            x = "sample",
            y = ylab,
            fill = paste(interestingGroups, collapse = ":\n"))

    if (is.numeric(warnLimit)) {
        p <- p + qcWarnLine(warnLimit)
    }

    if (is(fill, "ScaleDiscrete")) {
        p <- p + fill
    }

    if (isTRUE(flip)) {
        p <- p + coord_flip()
    }

    if (identical(interestingGroups, "sampleName")) {
        p <- p + guides(fill = FALSE)
    }

    p
}



.plot53Bias.bcbioRNASeq <- function(
    object,
    interestingGroups,
    warnLimit = 2L,
    fill = viridis::scale_fill_viridis(discrete = TRUE),
    flip = TRUE,
    title = TRUE) {
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }
    .plot53Bias(
        object = metrics(object),
        interestingGroups = interestingGroups,
        warnLimit = warnLimit,
        fill = fill,
        flip = flip,
        title = title)
}



# Methods ======================================================================
#' @rdname plot53Bias
#' @export
setMethod(
    "plot53Bias",
    signature("bcbioRNASeq"),
    .plot53Bias.bcbioRNASeq)



#' @rdname plot53Bias
#' @export
setMethod(
    "plot53Bias",
    signature("data.frame"),
    .plot53Bias)



# Aliases ======================================================================
#' @rdname plot53Bias
#' @export
plot5x3Bias <- function(...) {
    plot53Bias(...)
}
