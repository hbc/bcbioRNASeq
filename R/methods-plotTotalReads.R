#' Plot Total Reads
#'
#' @rdname plotTotalReads
#' @name plotTotalReads
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inheritParams general
#'
#' @param interestingGroups Category to use to group samples. In the plotting
#'   functions, this will define color and shape, where applicable. If unset,
#'   this is automatically determined by the metadata set inside the
#'   [bcbioRNASeq] object. When set to `NULL`, this will default to
#'   `sampleName`.
#' @param passLimit Threshold to plot pass color marker.
#' @param warnLimit Threshold to plot warning color marker.
#' @param fill Desired ggplot fill scale. Defaults to
#'   [scale_fill_viridis()]. Must supply discrete values. When set to
#'   `NULL`, the default ggplot2 color palette will be used. If manual color
#'   definitions are desired, we recommend using [scale_fill_manual()].
#' @param flip Flip x and y axes.
#'
#' @return [ggplot].
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq
#' plotTotalReads(bcb)
#' plotTotalReads(
#'     bcb,
#'     interestingGroups = "sampleName",
#'     fill = NULL)
#'
#' # data.frame
#' df <- metrics(bcb)
#' plotTotalReads(df)
NULL



# Constructors =================================================================
#' @importFrom bcbioBase uniteInterestingGroups
#' @importFrom ggplot2 aes_ coord_flip geom_bar ggplot guides labs
#' @importFrom rlang !!! syms
#' @importFrom tidyr unite
.plotTotalReads <- function(
    object,
    interestingGroups = "sampleName",
    passLimit = 20L,
    warnLimit = 10L,
    fill = scale_fill_viridis(discrete = TRUE),
    flip = TRUE,
    title = TRUE) {
    assert_is_data.frame(object)
    assert_formal_interesting_groups(object, interestingGroups)
    assert_is_an_implicit_integer(passLimit)
    assert_all_are_non_negative(passLimit)
    assert_is_an_implicit_integer(warnLimit)
    assert_all_are_non_negative(warnLimit)
    .assert_formal_scale_discrete(fill)
    assert_is_a_bool(flip)
    .assert_formal_title(title)

    if (isTRUE(title)) {
        title <- "total reads"
    } else if (!is.character(title)) {
        title <- NULL
    }

    data <- uniteInterestingGroups(object, interestingGroups)

    p <- ggplot(
        data = data,
        mapping = aes_(
            x = ~sampleName,
            y = ~totalReads / 1e6L,
            fill = ~interestingGroups)
    ) +
        geom_bar(stat = "identity") +
        labs(
            title = title,
            x = "sample",
            y = "total reads (million)",
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



.plotTotalReads.bcbioRNASeq <- function(  # nolint
    object,
    interestingGroups,
    passLimit = 20L,
    warnLimit = 10L,
    fill = scale_fill_viridis(discrete = TRUE),
    flip = TRUE,
    title = TRUE) {
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }
    .plotTotalReads(
        object = metrics(object),
        interestingGroups = interestingGroups,
        passLimit = passLimit,
        warnLimit = warnLimit,
        fill = fill,
        flip = flip,
        title = title)
}



# Methods ======================================================================
#' @rdname plotTotalReads
#' @export
setMethod(
    "plotTotalReads",
    signature("bcbioRNASeq"),
    .plotTotalReads.bcbioRNASeq)



#' @rdname plotTotalReads
#' @export
setMethod(
    "plotTotalReads",
    signature("data.frame"),
    .plotTotalReads)
