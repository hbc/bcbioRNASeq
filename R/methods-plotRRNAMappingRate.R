#' Plot Ribosomal RNA (rRNA) Mapping Rate
#'
#' @name plotRRNAMappingRate
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inherit plotTotalReads
#'
#' @examples
#' plotRRNAMappingRate(bcb_small)
NULL



# Constructors =================================================================
#' @importFrom bcbioBase uniteInterestingGroups
#' @importFrom ggplot2 aes_ coord_flip geom_bar ggplot guides labs
.plotRRNAMappingRate <- function(
    object,
    interestingGroups = "sampleName",
    warnLimit = 10L,
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
        title <- "rRNA mapping rate"
    } else if (!is_a_string(title)) {
        title <- NULL
    }

    # Fix for legacy camel variant mismatch (e.g. rRnaRate).
    if (!"rrnaRate" %in% colnames(object)) {
        warn(paste(
            "`rrnaRate` is missing from `metrics()` slot.",
            updateMsg
        ))
        col <- grep(
            pattern = "rrnarate",
            x = colnames(object),
            ignore.case = TRUE,
            value = TRUE
        )
        assert_is_a_string(col)
        object[["rrnaRate"]] <- object[[col]]
        object[[col]] <- NULL
    }

    data <- uniteInterestingGroups(object, interestingGroups)

    p <- ggplot(
        data = data,
        mapping = aes_(
            x = ~sampleName,
            y = ~rrnaRate * 100L,
            fill = ~interestingGroups
        )
    ) +
        geom_bar(stat = "identity") +
        labs(
            title = title,
            x = "sample",
            y = "rRNA mapping rate (%)",
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
#' @rdname plotRRNAMappingRate
#' @export
setMethod(
    "plotRRNAMappingRate",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        warnLimit = 10L,
        fill = scale_fill_viridis(discrete = TRUE),
        flip = TRUE,
        title = TRUE
    ) {
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        .plotRRNAMappingRate(
            object = metrics(object),
            interestingGroups = interestingGroups,
            warnLimit = warnLimit,
            fill = fill,
            flip = flip,
            title = title
        )
    }
)
