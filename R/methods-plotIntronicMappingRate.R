#' Plot Intronic Mapping Rate
#'
#' @name plotIntronicMappingRate
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inherit plotTotalReads
#'
#' @examples
#' plotIntronicMappingRate(bcb_small)
NULL



# Constructors =================================================================
#' @importFrom bcbioBase uniteInterestingGroups
#' @importFrom ggplot2 aes_ coord_flip geom_bar ggplot guides labs ylim
.plotIntronicMappingRate <- function(
    object,
    interestingGroups = "sampleName",
    warnLimit = 20L,
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
    if (isTRUE(title)) {
        title <- "intronic mapping rate"
    } else if (!is_a_string(title)) {
        title <- NULL
    }

    data <- uniteInterestingGroups(object, interestingGroups)

    p <- ggplot(
        data = data,
        mapping = aes_(
            x = ~sampleName,
            y = ~intronicRate * 100L,
            fill = ~interestingGroups
        )
    ) +
        geom_bar(stat = "identity") +
        labs(
            title = title,
            x = "sample",
            y = "intronic mapping rate (%)",
            fill = paste(interestingGroups, collapse = ":\n")
        ) +
        ylim(0L, 100L)

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
#' @rdname plotIntronicMappingRate
#' @export
setMethod(
    "plotIntronicMappingRate",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        warnLimit = 20L,
        fill = scale_fill_viridis(discrete = TRUE),
        flip = TRUE,
        title = TRUE
    ) {
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        .plotIntronicMappingRate(
            object = metrics(object),
            interestingGroups = interestingGroups,
            warnLimit = warnLimit,
            fill = fill,
            flip = flip,
            title = title
        )
    }
)
