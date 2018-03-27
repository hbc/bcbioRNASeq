#' Plot Exonic Mapping Rate
#'
#' @name plotExonicMappingRate
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @examples
#' plotExonicMappingRate(bcb_small)
NULL



# Constructors =================================================================
.plotExonicMappingRate <- function(
    object,
    interestingGroups,
    passLimit = 60L,
    fill = scale_fill_viridis(discrete = TRUE),
    flip = TRUE,
    title = "exonic mapping rate"
) {
    validObject(object)
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }
    assertIsAnImplicitInteger(passLimit)
    assert_all_are_non_negative(passLimit)
    assertIsFillScaleDiscreteOrNULL(fill)
    assert_is_a_bool(flip)
    assertIsAStringOrNULL(title)

    metrics <- metrics(object) %>%
        uniteInterestingGroups(interestingGroups)

    p <- ggplot(
        data = metrics,
        mapping = aes_(
            x = ~sampleName,
            y = ~exonicRate * 100L,
            fill = ~interestingGroups
        )
    ) +
        geom_bar(stat = "identity") +
        labs(
            title = title,
            x = "sample",
            y = "exonic mapping rate (%)",
            fill = paste(interestingGroups, collapse = ":\n")
        ) +
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
    .plotExonicMappingRate
)
