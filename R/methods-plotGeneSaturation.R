#' Plot Gene Detection Saturation
#'
#' @rdname plotGeneSaturation
#' @name plotGeneSaturation
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inherit plotTotalReads
#' @inheritParams plotCountDensity
#' @inheritParams plotGene
#'
#' @param counts Object containing a count matrix.
#' @param minCounts Numeric value for filtering the counts matrix before
#'   plotting.
#' @param trendline Include a trendline for each group.
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq
#' plotGeneSaturation(bcb)
#' plotGeneSaturation(
#'     bcb,
#'     interestingGroups = "sampleName",
#'     color = NULL)
#'
#' # data.frame, matrix
#' df <- metrics(bcb)
#' counts <- counts(bcb, normalized = "rlog")
#' plotGeneSaturation(df, counts = counts)
NULL



# Constructors =================================================================
#' @importFrom bcbioBase uniteInterestingGroups
#' @importFrom ggplot2 aes_ geom_point geom_smooth ggplot labs
.plotGeneSaturation <- function(
    object,
    counts,
    interestingGroups = "sampleName",
    minCounts = 0L,
    trendline = TRUE,
    color = scale_color_viridis(discrete = TRUE),
    title = TRUE) {
    assert_is_data.frame(object)
    assert_formal_interesting_groups(object, interestingGroups)
    assert_is_an_implicit_integer(minCounts)
    assert_all_are_non_negative(minCounts)
    assert_is_a_bool(trendline)
    .assert_formal_scale_discrete(color)
    .assert_formal_title(title)

    if (isTRUE(title)) {
        title <- "gene saturation"
    } else if (!is.character(title)) {
        title <- NULL
    }

    data <- uniteInterestingGroups(object, interestingGroups)

    p <- ggplot(
        data = data,
        mapping = aes_(
            x = ~mappedReads / 1e6L,
            y = colSums(counts > minCounts),
            color = ~interestingGroups)
    ) +
        geom_point(size = 3L) +
        labs(
            title = title,
            x = "mapped reads (million)",
            y = "genes",
            color = paste(interestingGroups, collapse = ":\n")
        )

    if (isTRUE(trendline)) {
        p <- p + geom_smooth(method = "lm", se = FALSE)
    }

    if (is(color, "ScaleDiscrete")) {
        p <- p + color
    }

    p
}



.plotGeneSaturation.bcbioRNASeq <- function(  # nolint
    object,
    interestingGroups,
    normalized = "tmm",
    minCounts = 0L,
    trendline = TRUE,
    color = scale_color_viridis(discrete = TRUE),
    title = TRUE) {
    assert_is_a_string(normalized)
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }
    .plotGeneSaturation(
        object = metrics(object),
        counts = counts(object, normalized = normalized),
        interestingGroups = interestingGroups,
        minCounts = minCounts,
        trendline = trendline,
        color = color,
        title = title)
}



# Methods ======================================================================
#' @rdname plotGeneSaturation
#' @export
setMethod(
    "plotGeneSaturation",
    signature(
        object = "bcbioRNASeq",
        counts = "missing"),
    .plotGeneSaturation.bcbioRNASeq)



#' @rdname plotGeneSaturation
#' @export
setMethod(
    "plotGeneSaturation",
    signature(
        object = "data.frame",
        counts = "matrix"),
    .plotGeneSaturation)
