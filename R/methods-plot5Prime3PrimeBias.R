#' Plot 5'->3' Bias
#'
#' @name plot5Prime3PrimeBias
#' @family Quality Control Plots
#' @author Michael Steinbaugh
#'
#' @inherit plotTotalReads
#'
#' @examples
#' plot5Prime3PrimeBias(bcb_small)
NULL



# Constructors =================================================================
#' @importFrom bcbioBase uniteInterestingGroups
#' @importFrom ggplot2 aes_string coord_flip geom_bar ggplot guides labs
.plot5Prime3PrimeBias <- function(
    object,
    interestingGroups = "sampleName",
    warnLimit = 2L,
    fill = scale_fill_viridis(discrete = TRUE),
    flip = TRUE,
    title = TRUE
) {
    assert_is_data.frame(object)
    assertFormalInterestingGroups(object, interestingGroups)
    assertIsAnImplicitInteger(warnLimit)
    assert_all_are_non_negative(warnLimit)
    assertIsFillScaleDiscreteOrNULL(fill)
    assert_is_a_bool(flip)

    # Title
    ylab <- "5'->3' bias"
    if (isTRUE(title)) {
        title <- ylab
    } else if (!is_a_string(title)) {
        title <- NULL
    }

    data <- uniteInterestingGroups(object, interestingGroups)

    # Legacy code: make sure `x53Bias` is camel sanitized to `x5x3Bias`.
    # The internal camel method has been updated in basejump 0.1.11.
    if ("x53Bias" %in% colnames(data)) {
        data <- dplyr::rename(data, "x5x3Bias" = "x53Bias")
    }

    p <- ggplot(
        data = data,
        mapping = aes_string(
            x = "sampleName",
            y = "x5x3Bias",
            fill = "interestingGroups"
        )
    ) +
        geom_bar(stat = "identity") +
        labs(
            title = title,
            x = "sample",
            y = ylab,
            fill = paste(interestingGroups, collapse = ":\n")
        )

    if (is_positive(warnLimit)) {
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



# Methods ======================================================================
#' @rdname plot5Prime3PrimeBias
#' @export
setMethod(
    "plot5Prime3PrimeBias",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        warnLimit = 2L,
        fill = scale_fill_viridis(discrete = TRUE),
        flip = TRUE,
        title = TRUE) {
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        .plot5Prime3PrimeBias(
            object = metrics(object),
            interestingGroups = interestingGroups,
            warnLimit = warnLimit,
            fill = fill,
            flip = flip,
            title = title
        )
    }
)
