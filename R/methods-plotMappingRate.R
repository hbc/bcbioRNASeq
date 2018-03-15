#' Plot Mapping Rate
#'
#' @name plotMappingRate
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inherit plotTotalReads
#'
#' @examples
#' plotMappingRate(bcb_small)
NULL



# Constructors =================================================================
#' @importFrom bcbioBase uniteInterestingGroups
#' @importFrom ggplot2 aes_ coord_flip geom_bar ggplot guides labs ylim
.plotMappingRate <- function(
    object,
    interestingGroups = "sampleName",
    passLimit = 90L,
    warnLimit = 70L,
    fill = scale_fill_viridis(discrete = TRUE),
    flip = TRUE,
    title = TRUE
) {
    assert_is_data.frame(object)
    assertFormalInterestingGroups(object, interestingGroups)
    assertIsAnImplicitInteger(passLimit)
    assert_all_are_non_negative(passLimit)
    assertIsAnImplicitInteger(warnLimit)
    assert_all_are_non_negative(warnLimit)
    assertIsFillScaleDiscreteOrNULL(fill)
    assert_is_a_bool(flip)

    # Title
    if (isTRUE(title)) {
        title <- "mapping rate"
    } else if (!is_a_string(title)) {
        title <- NULL
    }

    data <- uniteInterestingGroups(object, interestingGroups)

    p <- ggplot(
        data = data,
        mapping = aes_(
            x = ~sampleName,
            y = ~mappedReads / totalReads * 100L,
            fill = ~interestingGroups
        )
    ) +
        geom_bar(stat = "identity") +
        ylim(0L, 100L) +
        labs(
            title = title,
            x = "sample",
            y = "mapping rate (%)",
            fill = paste(interestingGroups, collapse = ":\n")
        )

    if (is_positive(passLimit)) {
        p <- p + qcPassLine(passLimit)
    }
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
#' @rdname plotMappingRate
#' @export
setMethod(
    "plotMappingRate",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        passLimit = 90L,
        warnLimit = 70L,
        fill = scale_fill_viridis(discrete = TRUE),
        flip = TRUE,
        title = TRUE
    ) {
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        .plotMappingRate(
            object = metrics(object),
            interestingGroups = interestingGroups,
            passLimit = passLimit,
            warnLimit = warnLimit,
            fill = fill,
            flip = flip,
            title = title
        )
    }
)
