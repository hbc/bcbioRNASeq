#' Plot Volcano
#'
#' @name plotVolcano
#' @family Differential Expression Functions
#' @author John Hutchinson, Michael Steinbaugh, Lorena Pantano
#'
#' @inheritParams bcbioBase::plotHeatmap
#' @inheritParams general
#' @param direction Plot "`both`", "`up`", or "`down`" directions.
#' @param histograms Show LFC and P value histograms.
#'
#' @seealso This function is an updated variant of
#'   `CHBUtils::volcano_density_plot()`.
#'
#' @return `ggplot`.
#'
#' @examples
#' gene2symbol <- gene2symbol(bcb_small)
#'
#' # DESeqResults ====
#' # Color DEGs in each direction separately
#' plotVolcano(
#'     object = res_small,
#'     sigPointColor = c(
#'         upregulated = "purple",
#'         downregulated = "orange"
#'     )
#' )
#'
#' # Label DEGs with a single color
#' plotVolcano(res_small, sigPointColor = "purple")
#'
#' # Directional support
#' plotVolcano(
#'     object = res_small,
#'     direction = "up",
#'     ntop = 5L,
#'     gene2symbol = gene2symbol,
#'     histograms = TRUE
#' )
#' plotVolcano(
#'     object = res_small,
#'     direction = "down",
#'     ntop = 5L,
#'     gene2symbol = gene2symbol,
#'     histograms = TRUE
#' )
#'
#' # Return coordinates as a data.frame
#' x <- plotVolcano(res_small, return = "data.frame")
#' glimpse(x)
NULL



# Methods ======================================================================
#' @rdname plotVolcano
#' @export
setMethod(
    "plotVolcano",
    signature("DESeqResults"),
    function(
        object,
        alpha,
        lfc = 0L,
        genes = NULL,
        gene2symbol = NULL,
        ntop = 0L,
        direction = c("both", "up", "down"),
        pointColor = "gray50",
        sigPointColor = c(
            upregulated = "purple",
            downregulated = "orange"
        ),
        histograms = FALSE,
        return = c("ggplot", "data.frame")
    ) {
        validObject(object)
        if (missing(alpha)) {
            alpha <- metadata(object)[["alpha"]]
        }
        assert_all_are_in_left_open_range(
            x = alpha,
            lower = 0L,
            upper = 1L
        )
        assert_is_a_number(lfc)
        assertFormalGene2symbol(object, genes, gene2symbol)
        assertIsImplicitInteger(ntop)
        direction <- match.arg(direction)
        assert_all_are_non_negative(c(lfc, ntop))
        assert_is_a_string(pointColor)
        assert_is_character(sigPointColor)
        if (is_a_string(sigPointColor)) {
            sigPointColor <- c(
                "upregulated" = sigPointColor,
                "downregulated" = sigPointColor
            )
        }
        assert_is_of_length(sigPointColor, 2L)
        assert_is_a_bool(histograms)
        return <- match.arg(return)

        data <- object %>%
            as.data.frame() %>%
            rownames_to_column("geneID") %>%
            as_tibble() %>%
            camel() %>%
            # Select columns used for plots
            .[, c("geneID", "baseMean", "log2FoldChange", "padj")] %>%
            # Remove rows with any NA values (e.g. padj)
            .[complete.cases(.), ] %>%
            # Remove rows with zero counts
            .[.[["baseMean"]] > 0L, , drop = FALSE] %>%
            # Negative log10 transform the P values. Add `1e-10` here to prevent
            # `Inf` values resulting from log transformation. Consequently, this
            # will gate the upper bound of the y-axis at 10.
            mutate(negLog10Pvalue = -log10(!!sym("padj") + 1e-10)) %>%
            # Calculate rank score. This is used to determine the `ntop` order.
            mutate(
                rankScore = !!sym("negLog10Pvalue") *
                    abs(!!sym("log2FoldChange"))
            ) %>%
            arrange(desc(!!sym("rankScore"))) %>%
            mutate(rank = row_number()) %>%
            .degColors(alpha = alpha, lfcThreshold = lfc)

        if (direction == "up") {
            data <- data[data[["log2FoldChange"]] > 0L, ]
        } else if (direction == "down") {
            data <- data[data[["log2FoldChange"]] < 0L, ]
        }

        # Gene-to-symbol mappings
        if (is.data.frame(gene2symbol)) {
            assertIsGene2symbol(gene2symbol)
            labelCol <- "geneName"
            data <- left_join(data, gene2symbol, by = "geneID")
        } else {
            labelCol <- "geneID"
        }

        # Early return data frame, if desired
        if (return == "data.frame") {
            return(data)
        }

        # LFC density ==========================================================
        lfcHist <- ggplot(
            data,
            mapping = aes_string(x = "log2FoldChange")
        ) +
            geom_density(
                color = NA,
                fill = pointColor
            ) +
            labs(x = "log2 fold change") +
            theme(
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.title.y = element_blank(),
                legend.position = "none"
            )

        # P value density ======================================================
        pvalueHist <- ggplot(
            data = data,
            mapping = aes_string(x = "negLog10Pvalue")
        ) +
            geom_density(
                color = NA,
                fill = pointColor
            ) +
            labs(x = "-log10 adj p value") +
            theme(
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.title.y = element_blank(),
                legend.position = "none"
            )

        # Volcano plot =========================================================
        p <- ggplot(
            data = data,
            mapping = aes_string(
                x = "log2FoldChange",
                y = "negLog10Pvalue",
                color = "color"
            )
        ) +
            geom_vline(
                xintercept = 0L,
                size = 0.5,
                color = pointColor
            ) +
            geom_point() +
            guides(color = FALSE) +
            labs(
                title = contrastName(object),
                x = "log2 fold change",
                y = "-log10 adj p value"
            )

        if (is_a_string(pointColor) && is.character(sigPointColor)) {
            p <- p +
                scale_color_manual(
                    values = c(
                        "nonsignificant" = pointColor,
                        "upregulated" = sigPointColor[[1L]],
                        "downregulated" = sigPointColor[[2L]]
                    )
                )
        }

        # Gene text labels =====================================================
        if (is.null(genes) && is_positive(ntop)) {
            genes <- data[1L:ntop, "geneID", drop = TRUE]
        }
        if (is.character(genes)) {
            assert_is_subset(genes, data[["geneID"]])
            labelData <- data[data[["geneID"]] %in% genes, ]
            p <- p +
                .geomLabel(
                    data = labelData,
                    mapping = aes_string(
                        x = "log2FoldChange",
                        y = "negLog10Pvalue",
                        label = labelCol
                    )
                )
        }

        # Grid layout ==========================================================
        if (isTRUE(histograms)) {
            ggdraw() +
                # Coordinates are relative to lower left corner
                draw_plot(
                    plot = p,
                    x = 0L, y = 0.15,
                    width = 1L, height = 0.85
                ) +
                draw_plot(
                    plot = lfcHist,
                    x = 0L, y = 0L,
                    width = 0.5, height = 0.15
                ) +
                draw_plot(
                    plot = pvalueHist,
                    x = 0.5, y = 0L,
                    width = 0.5, height = 0.15
                )
        } else {
            p
        }
    }
)
