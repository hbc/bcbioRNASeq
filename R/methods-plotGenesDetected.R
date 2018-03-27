#' Plot Genes Detected
#'
#' @name plotGenesDetected
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @examples
#' # Minimal example with fewer genes, so disable limits
#' plotGenesDetected(bcb_small, passLimit = 0L, warnLimit = 0L)
NULL



# Constructors =================================================================
.plotGenesDetected <- function(
    object,
    interestingGroups,
    passLimit = 20000L,
    warnLimit = 15000L,
    minCounts = 0L,
    fill = scale_fill_viridis(discrete = TRUE),
    flip = TRUE,
    title = "genes detected"
) {
    validObject(object)
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }
    assertIsAnImplicitInteger(passLimit)
    assert_all_are_non_negative(passLimit)
    assertIsAnImplicitInteger(warnLimit)
    assert_all_are_non_negative(warnLimit)
    assertIsAnImplicitInteger(minCounts)
    assert_all_are_non_negative(minCounts)
    assertIsFillScaleDiscreteOrNULL(fill)
    assert_is_a_bool(flip)
    assertIsAStringOrNULL(title)

    metrics <- metrics(object) %>%
        uniteInterestingGroups(interestingGroups)
    counts <- counts(object, normalized = FALSE)

    p <- ggplot(
        data = metrics,
        mapping = aes_(
            x = ~sampleName,
            y = colSums(counts > minCounts),
            fill = ~interestingGroups
        )
    ) +
        geom_bar(stat = "identity") +
        labs(
            title = title,
            x = "sample",
            y = "gene count",
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
#' @rdname plotGenesDetected
#' @export
setMethod(
    "plotGenesDetected",
    signature("bcbioRNASeq"),
    .plotGenesDetected
)
