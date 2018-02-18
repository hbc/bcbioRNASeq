#' Plot Intronic Mapping Rate
#'
#' @rdname plotIntronicMappingRate
#' @name plotIntronicMappingRate
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inherit plotTotalReads
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq
#' plotIntronicMappingRate(bcb)
#' plotIntronicMappingRate(
#'     bcb,
#'     interestingGroups = "sampleName",
#'     fill = NULL)
#'
#' # data.frame
#' df <- metrics(bcb)
#' plotIntronicMappingRate(df)
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
    title = TRUE) {
    assert_is_data.frame(object)
    assert_formal_interesting_groups(object, interestingGroups)
    assert_is_an_implicit_integer(warnLimit)
    assert_all_are_non_negative(warnLimit)
    .assert_formal_scale_discrete(fill)
    assert_is_a_bool(flip)
    .assert_formal_title(title)

    if (isTRUE(title)) {
        title <- "intronic mapping rate"
    } else if (!is.character(title)) {
        title <- NULL
    }

    data <- uniteInterestingGroups(object, interestingGroups)

    p <- ggplot(
        data = data,
        mapping = aes_(
            x = ~sampleName,
            y = ~intronicRate * 100L,
            fill = ~interestingGroups)
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



.plotIntronicMappingRate.bcbioRNASeq <- function(  # nolint
    object,
    interestingGroups,
    warnLimit = 20L,
    fill = scale_fill_viridis(discrete = TRUE),
    flip = TRUE,
    title = TRUE) {
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }
    .plotIntronicMappingRate(
        object = metrics(object),
        interestingGroups = interestingGroups,
        warnLimit = warnLimit,
        fill = fill,
        flip = flip,
        title = title)
}



# Methods ======================================================================
#' @rdname plotIntronicMappingRate
#' @export
setMethod(
    "plotIntronicMappingRate",
    signature("bcbioRNASeq"),
    .plotIntronicMappingRate.bcbioRNASeq)



#' @rdname plotIntronicMappingRate
#' @export
setMethod(
    "plotIntronicMappingRate",
    signature("data.frame"),
    .plotIntronicMappingRate)
