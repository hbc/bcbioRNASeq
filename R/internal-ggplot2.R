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



# ggproto objects ==============================================================
geneMedianLine <- stat_summary(
    fun.y = median,
    fun.ymin = median,
    fun.ymax = median,
    geom = "crossbar",
    show.legend = FALSE,
    width = 0.5
)



# geom functions ==============================================================
geomLabel <- function(data, mapping, size = 4L) {
    geom_label_repel(
        data = data,
        mapping = mapping,
        arrow = arrow(length = unit(0.01, "npc")),
        box.padding = unit(0.5, "lines"),
        fontface = "bold",
        force = 1L,
        point.padding = unit(0.75, "lines"),
        segment.size = 0.5,
        show.legend = FALSE,
        size = size
    )
}



genePoint <- function(size = 3L, alpha = 1L) {
    geom_point(
        size = size,
        alpha = alpha,
        position = position_jitterdodge(dodge.width = 0.9)
    )
}



qcPassLine <- function(intercept) {
    assert_is_a_number(intercept)
    assert_all_are_non_negative(intercept)
    geom_hline(
        alpha = qcLineAlpha,
        color = qcPassColor,
        linetype = qcLineType,
        size = qcLineSize,
        yintercept = intercept
    )
}



qcWarnLine <- function(intercept) {
    assert_is_a_number(intercept)
    assert_all_are_non_negative(intercept)
    geom_hline(
        alpha = qcLineAlpha,
        color = qcWarnColor,
        linetype = qcLineType,
        size = qcLineSize,
        yintercept = intercept
    )
}



# Calculate a numeric vector to define the colors
# -1: downregulated
#  0: not significant
#  1: upregulated
.degColors <- function(data, alpha, lfcThreshold = 0L) {
    color <- mapply(
        lfc = data[["log2FoldChange"]],
        padj = data[["padj"]],
        FUN = function(lfc, padj) {
            if (lfc > lfcThreshold & padj < alpha) {
                "upregulated"
            } else if (lfc < lfcThreshold & padj < alpha) {
                "downregulated"
            } else {
                "nonsignificant"
            }
        },
        SIMPLIFY = TRUE,
        USE.NAMES = FALSE
    )
    color <- as.factor(color)
    data[["color"]] <- color
    data
}
