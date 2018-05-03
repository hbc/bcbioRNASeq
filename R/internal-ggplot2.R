lineColor <- "black"

# Quality control plot colors (from inferno)
qcPassColor <- "#000004"
qcWarnColor <- "#BB3754"
qcFailColor <- "#FCFFA4"
qcCutoffColor <- "black"

# Quality control label appearance
qcLabelAlpha <- 0.75
qcLabelColor <- "white"
qcLabelFill <- "black"
qcLabelFontface <- "bold"
qcLabelPadding <- grid::unit(0.2, "lines")
qcLabelSize <- NA

# Plot label separator
labelSep <- " : "



# ggproto objects ==============================================================
.geneMedianLine <- stat_summary(
    fun.y = median,
    fun.ymin = median,
    fun.ymax = median,
    geom = "crossbar",
    show.legend = FALSE,
    width = 0.5
)



# geom functions ==============================================================
# Calculate a numeric vector to define the colors
# -1: downregulated
#  0: not significant
#  1: upregulated
.addIsDECol <- function(
    data,
    testCol = "padj",
    alpha,
    lfcCol = "log2FoldChange",
    lfcThreshold = 0L
) {
    # test: P value or S value
    test <- data[[testCol]]
    # lfc: log2 fold change cutoff
    lfc <- data[[lfcCol]]
    isDE <- mapply(
        test = test,
        lfc = lfc,
        FUN = function(test, lfc) {
            if (any(is.na(c(test, lfc)))) {
                # nonsignificant
                0L
            } else if (test < alpha & lfc > lfcThreshold) {
                # upregulated
                1L
            } else if (test < alpha & lfc < -lfcThreshold) {
                # downregulated
                -1L
            } else {
                0L
            }
        },
        SIMPLIFY = TRUE,
        USE.NAMES = FALSE
    )
    isDE <- as.factor(isDE)
    data[["isDE"]] <- isDE
    data
}



.geomLabel <- function(data = NULL, mapping = NULL, size = 4L) {
    geom_label_repel(
        data = data,
        mapping = mapping,
        arrow = arrow(length = unit(0.01, "npc")),
        box.padding = unit(0.5, "lines"),
        fill = "white",
        fontface = "bold",
        force = 1L,
        point.padding = unit(0.75, "lines"),
        segment.size = 0.5,
        show.legend = FALSE,
        size = size
    )
}



.genePoint <- function(size = 3L, alpha = 1L) {
    geom_point(
        size = size,
        alpha = alpha,
        position = position_jitterdodge(dodge.width = 0.9)
    )
}



.qcLine <- function(
    intercept,
    alpha = 0.75,
    color = "black",
    linetype = "dashed",
    size = 1L
) {
    assert_is_a_number(intercept)
    assert_all_are_non_negative(intercept)
    assert_is_a_number(alpha)
    assert_all_are_in_left_open_range(alpha, lower = 0L, upper = 1L)
    assert_is_a_string(color)
    assert_is_a_string(linetype)
    assert_is_a_number(size)
    assert_all_are_positive(size)
    geom_hline(
        alpha = alpha,
        color = color,
        linetype = linetype,
        size = size,
        yintercept = intercept
    )
}
