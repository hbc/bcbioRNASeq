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
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq
#' plotExonicMappingRate(bcb)
#' plotExonicMappingRate(
#'     bcb,
#'     interestingGroups = "sampleName",
#'     fill = NULL)
#'
#' # data.frame
#' df <- metrics(bcb)
#' plotExonicMappingRate(df)
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
    assert_formal_interesting_groups(object, interestingGroups)
    assert_is_an_implicit_integer(passLimit)
    .assert_formal_scale_discrete(fill)
    assert_is_a_bool(flip)
    .assert_formal_title(title)

    data <- uniteInterestingGroups(object, interestingGroups)

    if (isTRUE(title)) {
        title <- "exonic mapping rate"
    } else if (!is.character(title)) {
        title <- NULL
    }

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

    if (is.numeric(passLimit)) {
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



.plotExonicMappingRate.bcbioRNASeq <- function(  # nolint
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
}



# Methods ======================================================================
#' @rdname plotExonicMappingRate
#' @export
setMethod(
    "plotExonicMappingRate",
    signature("bcbioRNASeq"),
    .plotExonicMappingRate.bcbioRNASeq)



#' @rdname plotExonicMappingRate
#' @export
setMethod(
    "plotExonicMappingRate",
    signature("data.frame"),
    .plotExonicMappingRate)
