#' @importFrom ggplot2 aes_ aes_string annotation_logticks coord_fixed
#'   coord_flip element_blank element_text expand_limits facet_wrap geom_bar
#'   geom_boxplot geom_density geom_hline geom_jitter geom_line geom_point
#'   geom_polygon geom_ribbon geom_smooth geom_text ggplot ggtitle guides labs
#'   scale_color_manual scale_x_continuous scale_x_log10 scale_y_log10 theme
#'   xlab xlim ylab ylim
#' @importFrom viridis inferno

lineColor <- "black"

# Quality control plot colors
qcColors <- inferno(3)
qcPassColor <- qcColors[[1]]
qcWarnColor <- qcColors[[2]]
qcFailColor <- qcColors[[3]]
qcCutoffColor <- "black"

# Quality control label appearance
qcLabelAlpha <- 0.75
qcLabelColor <- "white"
qcLabelFill <- "black"
qcLabelFontface <- "bold"
qcLabelPadding <- unit(0.2, "lines")
qcLabelSize <- NA

# Quality control line appearance
qcLineAlpha <- 0.75
qcLineSize <- 1
qcLineType <- "dashed"

# Plot label separator
labelSep <- ": "

# Line functions
qcPassLine <- function(intercept) {
    geom_hline(
        alpha = qcLineAlpha,
        color = qcPassColor,
        linetype = qcLineType,
        size = qcLineSize,
        yintercept = intercept)
}

qcWarnLine <- function(intercept) {
    geom_hline(
        alpha = qcLineAlpha,
        color = qcWarnColor,
        linetype = qcLineType,
        size = qcLineSize,
        yintercept = intercept)
}
