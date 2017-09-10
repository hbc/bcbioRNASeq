lineColor <- "black"
qcColors <- inferno(3L)

# Quality control plot colors
qcCutoffColor <- "black"
qcPassColor <- qcColors[[1L]]
qcWarnColor <- qcColors[[2L]]
qcFailColor <- qcColors[[3L]]

# Quality control label appearance
qcLabelAlpha <- 0.75
qcLabelColor <- "white"
qcLabelFill <- "black"
qcLabelFontface <- "bold"
qcLabelPadding <- unit(0.2, "lines")
qcLabelSize <- NA

# Quality control line appearance
qcLineAlpha <- 0.75
qcLineSize <- 1L
qcLineType <- "dashed"



qcPassLine <- function(intercept) {
    geom_hline(alpha = qcLineAlpha,
               color = qcPassColor,
               linetype = qcLineType,
               size = qcLineSize,
               yintercept = intercept)
}



qcWarnLine <- function(intercept) {
    geom_hline(alpha = qcLineAlpha,
               color = qcWarnColor,
               linetype = qcLineType,
               size = qcLineSize,
               yintercept = intercept)
}
