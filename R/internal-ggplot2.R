#' @importFrom grid unit
#' @importFrom viridis inferno

lineColor <- "black"

# Quality control plot colors
qcColors <- viridis::inferno(3L)
qcPassColor <- qcColors[[1L]]
qcWarnColor <- qcColors[[2L]]
qcFailColor <- qcColors[[3L]]
qcCutoffColor <- "black"

# Quality control label appearance
qcLabelAlpha <- 0.75
qcLabelColor <- "white"
qcLabelFill <- "black"
qcLabelFontface <- "bold"
qcLabelPadding <- grid::unit(0.2, "lines")
qcLabelSize <- NA

# Quality control line appearance
qcLineAlpha <- 0.75
qcLineSize <- 1L
qcLineType <- "dashed"

# Plot label separator
labelSep <- ": "



# Line functions ===============================================================
#' @importFrom ggplot2 stat_summary
geneMedianLine <- stat_summary(
    fun.y = median,
    fun.ymin = median,
    fun.ymax = median,
    geom = "crossbar",
    show.legend = FALSE,
    width = 0.67)


#' @importFrom ggplot2 geom_hline
qcPassLine <- function(intercept) {
    assert_is_a_number(intercept)
    assert_all_are_non_negative(intercept)
    geom_hline(
        alpha = qcLineAlpha,
        color = qcPassColor,
        linetype = qcLineType,
        size = qcLineSize,
        yintercept = intercept)
}

#' @importFrom ggplot2 geom_hline
qcWarnLine <- function(intercept) {
    assert_is_a_number(intercept)
    assert_all_are_non_negative(intercept)
    geom_hline(
        alpha = qcLineAlpha,
        color = qcWarnColor,
        linetype = qcLineType,
        size = qcLineSize,
        yintercept = intercept)
}
