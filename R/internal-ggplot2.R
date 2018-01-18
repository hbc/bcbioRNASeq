#' @importFrom grid unit
#' @importFrom viridis inferno

lineColor <- "black"

# Quality control plot colors
qcColors <- viridis::inferno(3)
qcPassColor <- qcColors[[1]]
qcWarnColor <- qcColors[[2]]
qcFailColor <- qcColors[[3]]
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
qcLineSize <- 1
qcLineType <- "dashed"

# Plot label separator
labelSep <- ": "



# Line functions ===============================================================
#' @importFrom ggplot2 geom_hline
qcPassLine <- function(intercept) {
    geom_hline(
        alpha = qcLineAlpha,
        color = qcPassColor,
        linetype = qcLineType,
        size = qcLineSize,
        yintercept = intercept)
}

#' @importFrom ggplot2 geom_hline
qcWarnLine <- function(intercept) {
    geom_hline(
        alpha = qcLineAlpha,
        color = qcWarnColor,
        linetype = qcLineType,
        size = qcLineSize,
        yintercept = intercept)
}
