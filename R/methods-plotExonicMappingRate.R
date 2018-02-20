#' Plot Exonic Mapping Rate
#'
#' @rdname plotExonicMappingRate
#' @name plotExonicMappingRate
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inherit plotTotalReads
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))
#' plotExonicMappingRate(bcb)
NULL



# Constructors =================================================================
#' @importFrom bcbioBase uniteInterestingGroups
#' @importFrom ggplot2 aes_ coord_flip geom_bar ggplot guides labs ylim
.plotExonicMappingRate <- function(
    object,
    interestingGroups = "sampleName",
    passLimit = 60L,
    fill = scale_fill_viridis(discrete = TRUE),
    flip = TRUE,
    title = TRUE) {
    assert_is_data.frame(object)
    assertFormalInterestingGroups(object, interestingGroups)
    assertIsAnImplicitInteger(passLimit)
    assert_all_are_non_negative(passLimit)
    assertIsFillScaleDiscreteOrNULL(fill)
    assert_is_a_bool(flip)

    # Title
    if (isTRUE(title)) {
        title <- "exonic mapping rate"
    } else if (!is_a_string(title)) {
        title <- NULL
    }

    data <- uniteInterestingGroups(object, interestingGroups)

    p <- ggplot(
        data = data,
        mapping = aes_(
            x = ~sampleName,
            y = ~exonicRate * 100L,
            fill = ~interestingGroups)
    ) +
        geom_bar(stat = "identity") +
        labs(
            title = title,
            x = "sample",
            y = "exonic mapping rate (%)",
            fill = paste(interestingGroups, collapse = ":\n")) +
        ylim(0L, 100L)

    if (is_positive(passLimit)) {
        p <- p + qcPassLine(passLimit)
    }

    if (is(fill, "ScaleDiscrete")) {
        p <- p + scale_fill_viridis(discrete = TRUE)
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
#' @rdname plotExonicMappingRate
#' @export
setMethod(
    "plotExonicMappingRate",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        passLimit = 60L,
        fill = scale_fill_viridis(discrete = TRUE),
        flip = TRUE,
        title = TRUE) {
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        .plotExonicMappingRate(
            object = metrics(object),
            interestingGroups = interestingGroups,
            passLimit = passLimit,
            fill = fill,
            flip = flip,
            title = title)
    })
