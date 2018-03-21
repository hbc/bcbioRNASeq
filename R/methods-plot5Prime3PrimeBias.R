#' Plot 5'->3' Bias
#'
#' @name plot5Prime3PrimeBias
#' @family Quality Control Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @examples
#' plot5Prime3PrimeBias(bcb_small)
NULL



# Constructors =================================================================
#' @importFrom bcbioBase interestingGroups uniteInterestingGroups
#' @importFrom ggplot2 aes_string coord_flip geom_bar ggplot guides labs
.plot5Prime3PrimeBias <- function(
    object,
    interestingGroups,
    warnLimit = 2L,
    fill = scale_fill_viridis(discrete = TRUE),
    flip = TRUE,
    title = "5'->3' bias"
) {
    validObject(object)
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }
    assertIsAnImplicitInteger(warnLimit)
    assert_all_are_non_negative(warnLimit)
    assertIsFillScaleDiscreteOrNULL(fill)
    assert_is_a_bool(flip)
    assertIsAStringOrNULL(title)

    metrics <- metrics(object) %>%
        uniteInterestingGroups(interestingGroups)

    # Legacy code: make sure `x53Bias` is camel sanitized to `x5x3Bias`.
    # The internal camel method has been updated in basejump 0.1.11.
    if ("x53Bias" %in% colnames(metrics)) {
        metrics[["x5x3Bias"]] <- metrics[["x53Bias"]]
    }

    p <- ggplot(
        data = metrics,
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
            y = "5'->3' bias",
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
    .plot5Prime3PrimeBias
)
