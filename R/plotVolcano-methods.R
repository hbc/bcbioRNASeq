#' @name plotVolcano
#' @inherit bioverbs::plotVolcano
#' @family Differential Expression Functions
#' @author John Hutchinson, Michael Steinbaugh, Lorena Pantano
#'
#' @inheritParams general
#' @param ylim `scalar numeric`. Upper boundary limit for y-axis. Helps preserve
#'   dynamic range for gene sets containing highly significant P values (e.g.
#'   `1e-100`).
#' @param histograms `boolean`. Show LFC and P value histograms.
#' @param ... Additional arguments.
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



#' @rdname plotVolcano
#' @name plotVolcano
#' @importFrom bioverbs plotVolcano
#' @usage plotVolcano(object, ...)
#' @export
NULL



plotVolcano.bcbioRNASeq <-  # nolint
    function(
        object,
        alpha,
        lfcThreshold = 0L,
        ylim = 1e-10,
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
        assert_is_a_number(lfcThreshold)
        assert_is_a_number(ylim)
        assert_all_are_in_range(
            x = ylim,
            lower = 1e-100,
            upper = 1e-3
        )
        assertFormalGene2symbol(object, genes, gene2symbol)
        assertIsImplicitInteger(ntop)
        direction <- match.arg(direction)
        assert_all_are_non_negative(c(lfcThreshold, ntop))
        assert_is_a_string(pointColor)
        assert_is_character(sigPointColor)
        if (is_a_string(sigPointColor)) {
            sigPointColor <- c(
                upregulated = sigPointColor,
                downregulated = sigPointColor
            )
        }
        assert_is_of_length(sigPointColor, 2L)
        assert_is_a_bool(histograms)
        return <- match.arg(return)

        # Check to see if we should use `sval` instead of `padj`
        if ("svalue" %in% names(object)) {
            testCol <- "svalue"  # nocov
        } else {
            testCol <- "padj"
        }

        lfcCol <- "log2FoldChange"
        negLogTestCol <- camel(paste("neg", "log10", testCol))

        data <- object %>%
            as.data.frame() %>%
            rownames_to_column("geneID") %>%
            as_tibble() %>%
            camel() %>%
            # Select columns used for plots
            .[, c("geneID", "baseMean", lfcCol, testCol)] %>%
            # Remove rows with any NA values (e.g. padj)
            .[complete.cases(.), , drop = FALSE] %>%
            # Remove rows with zero counts
            .[.[["baseMean"]] > 0L, , drop = FALSE] %>%
            # Negative log10 transform the test values. Add `ylim` here to
            # prevent `Inf` values resulting from log transformation.
            # This will also define the upper bound of the y-axis.
            # Then calculate the rank score, which is used for `ntop`.
            mutate(
                !!sym(negLogTestCol) := -log10(!!sym(testCol) + !!ylim),
                rankScore = !!sym(negLogTestCol) * abs(!!sym(lfcCol))
            ) %>%
            arrange(desc(!!sym("rankScore"))) %>%
            mutate(rank = row_number()) %>%
            .addIsDECol(
                testCol = testCol,
                alpha = alpha,
                lfcThreshold = lfcThreshold
            )

        if (direction == "up") {
            data <- data[data[[lfcCol]] > 0L, , drop = FALSE]
        } else if (direction == "down") {
            data <- data[data[[lfcCol]] < 0L, , drop = FALSE]
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
            data <- data %>%
                as.data.frame() %>%
                column_to_rownames("geneID")
            return(data)
        }

        # LFC density ----------------------------------------------------------
        lfcHist <- ggplot(
            data = data,
            mapping = aes(x = !!sym(lfcCol))
        ) +
            geom_density(
                color = NA,
                fill = pointColor
            ) +
            scale_x_continuous(
                breaks = pretty_breaks(),
                expand = c(0L, 0L)
            ) +
            scale_y_continuous(expand = c(0L, 0L)) +
            labs(
                x = "log2 fold change",
                y = NULL
            ) +
            guides(fill = FALSE) +
            theme(
                axis.line.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank()
            )

        # P value density ------------------------------------------------------
        pvalueHist <- ggplot(
            data = data,
            mapping = aes(x = !!sym(negLogTestCol))
        ) +
            geom_density(
                color = NA,
                fill = pointColor
            ) +
            scale_x_continuous(
                breaks = pretty_breaks(),
                expand = c(0L, 0L)
            ) +
            scale_y_continuous(expand = c(0L, 0L)) +
            labs(
                x = "-log10 adj p value",
                y = NULL
            ) +
            guides(fill = FALSE) +
            theme(
                axis.line.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank()
            )

        # Volcano plot ---------------------------------------------------------
        p <- ggplot(
            data = data,
            mapping = aes(
                x = !!sym(lfcCol),
                y = !!sym(negLogTestCol),
                color = !!sym("isDE")
            )
        ) +
            geom_vline(
                xintercept = 0L,
                size = 0.5,
                color = pointColor
            ) +
            geom_point() +
            scale_x_continuous(breaks = pretty_breaks()) +
            scale_y_continuous(breaks = pretty_breaks()) +
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
                        # nonsignificant
                        "0" = pointColor,
                        # upregulated
                        "1" = sigPointColor[[1L]],
                        # downregulated
                        "-1" = sigPointColor[[2L]]
                    )
                )
        }

        # Gene text labels -----------------------------------------------------
        if (is.null(genes) && is_positive(ntop)) {
            genes <- data[1L:ntop, "geneID", drop = TRUE]
        }
        if (is.character(genes)) {
            assert_is_subset(genes, data[["geneID"]])
            labelData <- data[data[["geneID"]] %in% genes, , drop = FALSE]
            p <- p +
                bcbio_geom_label_repel(
                    data = labelData,
                    mapping = aes(
                        x = !!sym(lfcCol),
                        y = !!sym(negLogTestCol),
                        label = !!sym(labelCol)
                    )
                )
        }

        # Grid layout ----------------------------------------------------------
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
                    width = 0.45, height = 0.15
                ) +
                draw_plot(
                    plot = pvalueHist,
                    x = 0.55, y = 0L,
                    width = 0.45, height = 0.15
                )
        } else {
            p
        }
    }



#' @rdname plotVolcano
#' @export
setMethod(
    f = "plotVolcano",
    signature = signature("DESeqResults"),
    definition = plotVolcano.bcbioRNASeq
)
