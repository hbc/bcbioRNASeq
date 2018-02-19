#' Plot Genes Detected
#'
#' @rdname plotGenesDetected
#' @name plotGenesDetected
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inherit plotTotalReads
#' @inheritParams plotGeneSaturation
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq
#' plotGenesDetected(
#'     bcb,
#'     passLimit = 0L,
#'     warnLimit = 0L)
#' plotGenesDetected(
#'     bcb,
#'     interestingGroups = "sampleName",
#'     fill = NULL,
#'     passLimit = 0L,
#'     warnLimit = 0L)
#'
#' # data.frame, DESeqDataSet
#' df <- metrics(bcb)
#' dds <- bcbio(bcb, "DESeqDataSet")
#' plotGenesDetected(df, counts = dds)
#'
#' # data.frame, matrix
#' counts <- counts(bcb, normalized = TRUE)
#' plotGenesDetected(df, counts = counts)
NULL



# Constructors =================================================================
#' @importFrom bcbioBase uniteInterestingGroups
#' @importFrom ggplot2 aes_ coord_flip geom_bar ggplot guides labs
.plotGenesDetected <- function(
    object,
    counts,
    interestingGroups = "sampleName",
    passLimit = 20000L,
    warnLimit = 15000L,
    minCounts = 0L,
    fill = scale_fill_viridis(discrete = TRUE),
    flip = TRUE,
    title = TRUE) {
    assert_is_data.frame(object)
    assert_is_matrix(counts)
    assertFormalIntersectingGroups(object, interestingGroups)
    assertIsAnImplicitInteger(passLimit)
    assert_all_are_non_negative(passLimit)
    assertIsAnImplicitInteger(warnLimit)
    assert_all_are_non_negative(warnLimit)
    assertIsAnImplicitInteger(minCounts)
    assert_all_are_non_negative(minCounts)
    assertIsScaleFillDiscreteOrNULL(fill)
    assert_is_a_bool(flip)
    
    # Title
    if (isTRUE(title)) {
        title <- "genes detected"
    } else if (!is_a_string(title)) {
        title <- NULL
    }

    data <- uniteInterestingGroups(object, interestingGroups)

    p <- ggplot(
        data = data,
        mapping = aes_(
            x = ~sampleName,
            y = colSums(counts > minCounts),
            fill = ~interestingGroups)
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
    signature(
        object = "bcbioRNASeq",
        counts = "missing"),
    function(
        object,
        interestingGroups,
        passLimit = 20000L,
        warnLimit = 15000L,
        minCounts = 0L,
        fill = scale_fill_viridis(discrete = TRUE),
        flip = TRUE) {
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        .plotGenesDetected(
            object = metrics(object),
            counts = assay(object),
            interestingGroups = interestingGroups,
            passLimit = passLimit,
            warnLimit = warnLimit,
            minCounts = minCounts,
            fill = fill,
            flip = flip)
    })



#' @rdname plotGenesDetected
#' @export
setMethod(
    "plotGenesDetected",
    signature(
        object = "data.frame",
        counts = "DESeqDataSet"),
    function(
        object,
        counts,
        interestingGroups = "sampleName",
        passLimit = 20000L,
        warnLimit = 15000L,
        minCounts = 0L,
        fill = scale_fill_viridis(discrete = TRUE),
        flip = TRUE,
        title = TRUE) {
        .plotGenesDetected(
            object,
            counts = assay(counts),
            interestingGroups = interestingGroups,
            passLimit = passLimit,
            warnLimit = warnLimit,
            minCounts = minCounts,
            fill = fill,
            flip = flip,
            title = title)
    })



#' @rdname plotGenesDetected
#' @export
setMethod(
    "plotGenesDetected",
    signature(
        object = "data.frame",
        counts = "matrix"),
    .plotGenesDetected)
