#' Plot Mapped Reads
#'
#' @name plotMappedReads
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inherit plotTotalReads
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))
#' plotMappedReads(bcb)
NULL



# Constructors =================================================================
#' @importFrom bcbioBase uniteInterestingGroups
#' @importFrom ggplot2 aes_ coord_flip geom_bar ggplot guides labs
.plotMappedReads <- function(
    object,
    interestingGroups = "sampleName",
    passLimit = 20L,
    warnLimit = 10L,
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
        title <- "mapped reads"
    } else if (!is_a_string(title)) {
        title <- NULL
    }

    data <- uniteInterestingGroups(object, interestingGroups)

    p <- ggplot(
        data = data,
        mapping = aes_(
            x = ~sampleName,
            y = ~mappedReads / 1e6L,
            fill = ~interestingGroups
        )
    ) +
        geom_bar(stat = "identity") +
        labs(
            title = title,
            x = "sample",
            y = "mapped reads (million)",
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
#' @rdname plotMappedReads
#' @export
setMethod(
    "plotMappedReads",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        passLimit = 20L,
        warnLimit = 10L,
        fill = scale_fill_viridis(discrete = TRUE),
        flip = TRUE,
        title = TRUE
    ) {
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        .plotMappedReads(
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
