#' Plot Volcano
#'
#' @rdname plotVolcano
#' @name plotVolcano
#'
#' @param alpha Alpha level cutoff used for coloring.
#' @param padj Use P values adjusted for multiple comparisions.
#' @param lfc Log fold change ratio (base 2) cutoff for coloring.
#' @param genes Character vector of gene symbols to label.
#' @param ntop Number of top genes to label.
#' @param direction Plot `up`, `down`, or `both` (**default**) directions.
#' @param shadeColor Shading color for bounding box.
#' @param shadeAlpha Shading transparency alpha.
#' @param pointColor Point color.
#' @param pointAlpha Point transparency alpha.
#' @param pointOutlineColor Point outline color.
#' @param histograms Show LFC and P value histograms.
#'
#' @seealso This function is an updated variant of
#'   `CHBUtils::volcano_density_plot()`.
#'
#' @return Volcano plot arranged as grid (`grid = TRUE`), or [show()]
#'   individual [ggplot] (`grid = FALSE`).
#'
#' @examples
#' data(res)
#' plotVolcano(res)
#' plotVolcano(res, genes = "Sulf1")
#' plotVolcano(res, padj = FALSE, alpha = 0.01, lfc = 4L)
#' plotVolcano(res, histograms = FALSE, ntop = 5L)
NULL



# Constructors ====
.plotVolcano <- function(
    object,
    alpha = 0.1,
    padj = TRUE,
    lfc = 1L,
    genes = NULL,
    ntop = 0L,
    direction = "both",
    shadeColor = "green",
    shadeAlpha = 0.25,
    pointColor = "gray",
    pointAlpha = 0.75,
    pointOutlineColor = "darkgray",
    histograms = TRUE) {
    df <- as.data.frame(object)

    if (!any(direction %in% c("both", "down", "up")) |
        length(direction) > 1L) {
        stop("Direction must be both, up, or down")
    }


    # Generate stats tibble ====
    stats <- df %>%
        rownames_to_column("ensgene") %>%
        as("tibble") %>%
        camel %>%
        .[, c("ensgene", "log2FoldChange", "pvalue", "padj")] %>%
        .[!is.na(.[["log2FoldChange"]]), ]
    g2s <- detectOrganism(stats[["ensgene"]][[1L]]) %>%
        annotable(format = "gene2symbol")
    stats <- left_join(stats, g2s, by = "ensgene")

    # Negative log10 transform the P values
    # Add `1e-10` here to prevent `Inf` values resulting from `log10()`
    if (isTRUE(padj)) {
        stats <- stats %>%
            mutate(negLog10Pvalue = -log10(.data[["padj"]] + 1e-10))
        pTitle <- "padj value"
    } else {
        stats <- stats %>%
            mutate(negLog10Pvalue = -log10(.data[["pvalue"]] + 1e-10))
        pTitle <- "p value"
    }

    stats <- stats %>%
        # Calculate rank score
        mutate(rankScore = .data[["negLog10Pvalue"]] *
                   abs(.data[["log2FoldChange"]])) %>%
        arrange(desc(!!sym("rankScore")))


    # Text labels ====
    if (!is.null(genes)) {
        volcanoText <- stats %>%
            .[.[["symbol"]] %in% genes, , drop = FALSE]
    } else if (ntop > 0L) {
        volcanoText <- stats[1L:ntop, , drop = FALSE]
    } else {
        volcanoText <- NULL
    }


    # Plot ranges ====
    # Get range of LFC and P values to set up plot borders
    rangeLFC <-
        c(floor(min(na.omit(stats[["log2FoldChange"]]))),
          ceiling(max(na.omit(stats[["log2FoldChange"]]))))
    rangeNegLog10Pvalue <-
        c(floor(min(na.omit(stats[["negLog10Pvalue"]]))),
          ceiling(max(na.omit(stats[["negLog10Pvalue"]]))))


    # LFC density histogram ====
    lfcDensity <- stats[["log2FoldChange"]] %>%
        na.omit %>%
        density
    lfcDensityDf <- data.frame(
        x = lfcDensity[["x"]],
        y = lfcDensity[["y"]])
    lfcHist <- stats %>%
        ggplot(aes_(x = ~log2FoldChange)) +
        geom_density() +
        scale_x_continuous(limits = rangeLFC) +
        labs(x = "log2 fold change",
             y = "") +
        # Don't label density y-axis
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
    if (direction == "both" | direction == "up") {
        lfcHist <- lfcHist +
            geom_ribbon(
                data = filter(lfcDensityDf, .data[["x"]] > !!lfc),
                aes_(x = ~x, ymax = ~y),
                ymin = 0L,
                fill = shadeColor,
                alpha = shadeAlpha)
    }
    if (direction == "both" | direction == "down") {
        lfcHist <- lfcHist +
            geom_ribbon(
                data = filter(lfcDensityDf, .data[["x"]] < -UQ(lfc)),
                aes_(x = ~x, ymax = ~y),
                ymin = 0L,
                fill = shadeColor,
                alpha = shadeAlpha)
    }


    # P value density plot ====
    pvalueDensity <- stats[["negLog10Pvalue"]] %>%
        na.omit %>%
        density
    pvalueDensityDf <-
        data.frame(x = pvalueDensity[["x"]],
                   y = pvalueDensity[["y"]])
    pvalueHist <- stats %>%
        ggplot(aes_(x = ~negLog10Pvalue)) +
        geom_density() +
        geom_ribbon(data = filter(pvalueDensityDf,
                                  .data[["x"]] > -log10(UQ(alpha) + 1e-10)),
                    aes_(x = ~x, ymax = ~y),
                    ymin = 0L,
                    fill = shadeColor,
                    alpha = shadeAlpha) +
        labs(x = paste("-log10", pTitle),
             y = "") +
        # Don't label density y-axis
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank())


    # Volcano plot ====
    volcano <- stats %>%
        ggplot(aes_(x = ~log2FoldChange,
                    y = ~negLog10Pvalue)) +
        labs(x = "log2 fold change",
             y = paste("-log10", pTitle)) +
        geom_point(
            alpha = pointAlpha,
            color = pointOutlineColor,
            fill = pointColor,
            pch = 21L) +
        theme(legend.position = "none") +
        scale_x_continuous(limits = rangeLFC)
    if (!is.null(volcanoText)) {
        volcano <- volcano +
            geom_text_repel(
                data = volcanoText,
                aes_(x = ~log2FoldChange,
                     y = ~negLog10Pvalue,
                     label = ~symbol),
                arrow = arrow(length = unit(0.01, "npc")),
                box.padding = unit(0.5, "lines"),
                color = "black",
                fontface = "bold",
                force = 1L,
                point.padding = unit(0.75, "lines"),
                segment.color = "gray",
                segment.size = 0.5,
                show.legend = FALSE,
                size = 4L)
    }
    if (direction == "both" | direction == "up") {
        volcanoPolyUp <- with(stats, data.frame(
            x = as.numeric(c(
                lfc,
                lfc,
                max(rangeLFC),
                max(rangeLFC))),
            y = as.numeric(c(
                -log10(alpha + 1e-10),
                max(rangeNegLog10Pvalue),
                max(rangeNegLog10Pvalue),
                -log10(alpha + 1e-10)))))
        volcano <- volcano +
            geom_polygon(
                data = volcanoPolyUp,
                aes_(x = ~x, y = ~y),
                fill = shadeColor,
                alpha = shadeAlpha)
    }
    if (direction == "both" | direction == "down") {
        volcanoPolyDown <- with(stats, data.frame(
            x = as.numeric(c(
                -lfc,
                -lfc,
                min(rangeLFC),
                min(rangeLFC))),
            y = as.numeric(c(
                -log10(alpha + 1e-10),
                max(rangeNegLog10Pvalue),
                max(rangeNegLog10Pvalue),
                -log10(alpha + 1e-10)))))
        volcano <- volcano +
            geom_polygon(
                data = volcanoPolyDown,
                aes_(x = ~x, y = ~y),
                fill = shadeColor,
                alpha = shadeAlpha)
    }


    # Grid layout ====
    if (isTRUE(histograms)) {
        ggdraw() +
            # Coordinates are relative to lower left corner
            draw_plot(
                lfcHist,
                x = 0L, y = 0.7, width = 0.5, height = 0.3) +
            draw_plot(
                pvalueHist,
                x = 0.5, y = 0.7, width = 0.5, height = 0.3) +
            draw_plot(
                volcano, x = 0L, y = 0L, width = 1L, height = 0.7)
    } else {
        volcano + ggtitle("volcano")
    }
}



# Methods ====
#' @rdname plotVolcano
#' @export
setMethod("plotVolcano", "DESeqResults", function(object, alpha = NULL, ...) {
    if (is.null(alpha)) {
        alpha <- metadata(object)[["alpha"]]
    }
    .plotVolcano(object, alpha = alpha, ...)
})
