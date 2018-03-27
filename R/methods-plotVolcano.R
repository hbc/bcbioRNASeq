# TODO Improve y-axis break consistency from 1-10

#' Plot Volcano
#'
#' @name plotVolcano
#' @family Differential Expression Functions
#' @author John Hutchinson, Michael Steinbaugh, Lorena Pantano
#'
#' @inheritParams plotHeatmap
#' @inheritParams general
#' @param padj Use P values adjusted for multiple comparisions.
#' @param ntop Number of top genes to label.
#' @param direction Plot "`both`", "`up`", or "`down`" directions.
#' @param pointColor Point color.
#' @param pointAlpha Point transparency alpha.
#' @param pointOutlineColor Point outline color.
#' @param shadeColor Shading color for bounding box.
#' @param shadeAlpha Shading transparency alpha.
#' @param labelColor Gene label color.
#' @param histograms Show LFC and P value histograms.
#'
#' @seealso This function is an updated variant of
#'   `CHBUtils::volcano_density_plot()`.
#'
#' @return `ggplot`.
#'
#' @examples
#' # DESeqResults ====
#' plotVolcano(
#'     object = res_small,
#'     ntop = 5L,
#'     gene2symbol = gene2symbol(bcb_small)
#' )
NULL



# Constructors =================================================================
.plotVolcano <- function(
    object,
    alpha,
    padj = TRUE,
    lfc = 0L,
    genes = NULL,
    gene2symbol = NULL,
    ntop = 0L,
    direction = c("both", "up", "down"),
    pointColor = "gray",
    pointAlpha = 0.75,
    pointOutlineColor = "darkgray",
    shadeColor = "green",
    shadeAlpha = 0.25,
    labelColor = "black",
    histograms = TRUE
) {
    validObject(object)
    if (missing(alpha)) {
        alpha <- metadata(object)[["alpha"]]
    }
    assert_is_a_number(alpha)
    assert_is_a_bool(padj)
    assert_is_a_number(lfc)
    assertFormalGene2symbol(object, genes, gene2symbol)
    assertIsImplicitInteger(ntop)
    direction <- match.arg(direction)
    assert_is_a_string(pointColor)
    assert_is_a_number(pointAlpha)
    assert_is_a_string(pointOutlineColor)
    assert_is_a_string(shadeColor)
    assert_is_a_number(shadeAlpha)
    assert_all_are_non_negative(c(lfc, ntop))
    assert_all_are_in_left_open_range(
        x = c(alpha, pointAlpha, shadeAlpha),
        lower = 0L, upper = 1L
    )
    assert_is_a_string(labelColor)
    assert_is_a_bool(histograms)

    # Generate data `tibble`
    data <- object %>%
        as.data.frame() %>%
        rownames_to_column("geneID") %>%
        as_tibble() %>%
        camel(strict = FALSE) %>%
        # Keep genes with non-zero counts
        .[.[["baseMean"]] > 0L, , drop = FALSE] %>%
        # Keep genes with a fold change
        .[!is.na(.[["log2FoldChange"]]), , drop = FALSE] %>%
        # Keep genes with a P value
        .[!is.na(.[["pvalue"]]), , drop = FALSE] %>%
        # Select columns used for plots
        .[, c("geneID", "log2FoldChange", "pvalue", "padj")]

    # Negative log10 transform the P values
    # Add `1e-10` here to prevent `Inf` values resulting from `log10()`
    if (isTRUE(padj)) {
        data <- data %>%
            # Keep genes with an adjusted P value
            .[!is.na(.[["padj"]]), , drop = FALSE] %>%
            # log10 transform
            mutate(negLog10Pvalue = -log10(!!sym("padj") + 1e-10))
        pvalTitle <- "adj p value"
    } else {
        data <- data %>%
            mutate(negLog10Pvalue = -log10(!!sym("pvalue") + 1e-10))
        pvalTitle <- "p value"
    }

    # Calculate rank score
    data <- data %>%
        mutate(
            rankScore = !!sym("negLog10Pvalue") * abs(!!sym("log2FoldChange"))
        ) %>%
        arrange(desc(!!sym("rankScore")))

    # Gene text labels =========================================================
    if (is.null(genes) && is_positive(ntop)) {
        genes <- data[1L:ntop, "geneID", drop = TRUE]
    }
    if (is.character(genes)) {
        assert_is_subset(genes, data[["geneID"]])
        volcanoText <- data %>%
            .[.[["geneID"]] %in% genes, , drop = FALSE]
        if (is.data.frame(gene2symbol)) {
            labelCol <- "geneName"
            assertIsGene2symbol(gene2symbol)
            volcanoText <- left_join(volcanoText, gene2symbol, by = "geneID")
        } else {
            labelCol <- "geneID"
        }
    } else {
        volcanoText <- NULL
    }

    # Plot ranges ==============================================================
    # Get range of LFC and P values to set up plot borders
    rangeLFC <- c(
        floor(min(na.omit(data[["log2FoldChange"]]))),
        ceiling(max(na.omit(data[["log2FoldChange"]])))
    )
    rangeNegLog10Pvalue <- c(
        floor(min(na.omit(data[["negLog10Pvalue"]]))),
        ceiling(max(na.omit(data[["negLog10Pvalue"]])))
    )

    # LFC density ==============================================================
    lfcDensity <- data[["log2FoldChange"]] %>%
        na.omit() %>%
        density()
    lfcDensityDf <- data.frame(
        x = lfcDensity[["x"]],
        y = lfcDensity[["y"]]
    )
    lfcHist <- ggplot(
        data,
        mapping = aes_string(x = "log2FoldChange")
    ) +
        geom_density() +
        scale_x_continuous(limits = rangeLFC) +
        labs(x = "log2 fold change") +
        theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
        )
    if (direction == "both" || direction == "up") {
        lfcHist <- lfcHist +
            geom_ribbon(
                data = lfcDensityDf %>%
                    .[.[["x"]] > lfc, , drop = FALSE],
                mapping = aes_string(x = "x", ymax = "y"),
                ymin = 0L,
                fill = shadeColor,
                alpha = shadeAlpha
            )
    }
    if (direction == "both" || direction == "down") {
        lfcHist <- lfcHist +
            geom_ribbon(
                data = lfcDensityDf %>%
                    .[.[["x"]] < -lfc, ],
                mapping = aes_string(x = "x", ymax = "y"),
                ymin = 0L,
                fill = shadeColor,
                alpha = shadeAlpha
            )
    }

    # P value density ==========================================================
    pvalueDensity <- data[["negLog10Pvalue"]] %>%
        na.omit() %>%
        density()
    pvalueDensityDf <- data.frame(
        x = pvalueDensity[["x"]],
        y = pvalueDensity[["y"]]
    )
    pvalueHist <- ggplot(
        data,
        mapping = aes_string(x = "negLog10Pvalue")
    ) +
        geom_density() +
        geom_ribbon(
            data = pvalueDensityDf %>%
                .[.[["x"]] > -log10(alpha + 1e-10), , drop = FALSE],
            mapping = aes_string(x = "x", ymax = "y"),
            ymin = 0L,
            fill = shadeColor,
            alpha = shadeAlpha
        ) +
        labs(x = paste("-log10", pvalTitle)) +
        theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
        )

    # Volcano plot =============================================================
    volcano <- ggplot(
        data,
        mapping = aes_string(
            x = "log2FoldChange",
            y = "negLog10Pvalue"
        )
    ) +
        labs(
            x = "log2 fold change",
            y = paste("-log10", pvalTitle)
        ) +
        geom_point(
            alpha = pointAlpha,
            color = pointOutlineColor,
            fill = pointColor,
            pch = 21L
        ) +
        theme(legend.position = "none") +
        scale_x_continuous(limits = rangeLFC)
    if (is.data.frame(volcanoText)) {
        volcano <- volcano +
            geom_text_repel(
                data = volcanoText,
                mapping = aes_string(
                    x = "log2FoldChange",
                    y = "negLog10Pvalue",
                    label = labelCol
                ),
                arrow = arrow(length = unit(0.01, "npc")),
                box.padding = unit(0.5, "lines"),
                color = labelColor,
                fontface = "bold",
                force = 1L,
                point.padding = unit(0.75, "lines"),
                segment.color = labelColor,
                segment.size = 0.5,
                show.legend = FALSE,
                size = 4L
            )
    }
    if (direction == "both" || direction == "up") {
        volcanoPolyUp <- with(
            data,
            expr = data.frame(
                x = as.numeric(c(
                    lfc,
                    lfc,
                    max(rangeLFC),
                    max(rangeLFC)
                )),
                y = as.numeric(c(
                    -log10(alpha + 1e-10),
                    max(rangeNegLog10Pvalue),
                    max(rangeNegLog10Pvalue),
                    -log10(alpha + 1e-10)
                ))
            )
        )
        volcano <- volcano +
            geom_polygon(
                data = volcanoPolyUp,
                mapping = aes_string(x = "x", y = "y"),
                fill = shadeColor,
                alpha = shadeAlpha
            )
    }
    if (direction == "both" || direction == "down") {
        volcanoPolyDown <- with(
            data,
            expr = data.frame(
                x = as.numeric(c(
                    -lfc,
                    -lfc,
                    min(rangeLFC),
                    min(rangeLFC)
                )),
                y = as.numeric(c(
                    -log10(alpha + 1e-10),
                    max(rangeNegLog10Pvalue),
                    max(rangeNegLog10Pvalue),
                    -log10(alpha + 1e-10)
                ))
            )
        )
        volcano <- volcano +
            geom_polygon(
                data = volcanoPolyDown,
                mapping = aes_string(
                    x = "x",
                    y = "y"
                ),
                fill = shadeColor,
                alpha = shadeAlpha
            )
    }

    # Grid layout ==============================================================
    if (isTRUE(histograms)) {
        ggdraw() +
            # Coordinates are relative to lower left corner
            draw_plot(
                plot = lfcHist,
                x = 0L, y = 0.7,
                width = 0.5, height = 0.3
            ) +
            draw_plot(
                plot = pvalueHist,
                x = 0.5, y = 0.7,
                width = 0.5, height = 0.3
            ) +
            draw_plot(
                plot = volcano,
                x = 0L, y = 0L,
                width = 1L, height = 0.7
            )
    } else {
        volcano + ggtitle("volcano")
    }
}



# Methods ======================================================================
#' @rdname plotVolcano
#' @export
setMethod(
    "plotVolcano",
    signature("DESeqResults"),
    .plotVolcano
)
