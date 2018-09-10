#' Plot Volcano
#'
#' @name plotVolcano
#' @family Differential Expression Functions
#' @author John Hutchinson, Michael Steinbaugh, Lorena Pantano
#'
#' @inheritParams general
#' @param ylim `scalar numeric`. Upper boundary limit for y-axis. Helps preserve
#'   dynamic range for gene sets containing highly significant P values (e.g.
#'   `1e-100`).
#' @param histograms `boolean`. Show LFC and P value histograms.
#'
#' @seealso This function is an updated variant of
#'   `CHBUtils::volcano_density_plot()`.
#'
#' @return `ggplot`.
#'
#' @examples
#' gene2symbol <- gene2symbol(bcb_small)
#' print(gene2symbol)
#'
#' geneIDs <- head(gene2symbol[["geneID"]])
#' print(geneIDs)
#'
#' geneNames <- head(gene2symbol[["geneName"]])
#' print(geneNames)
#'
#' # DESeqResults ====
#' summary(res_small)
#'
#' # Color DEGs in each direction separately.
#' plotVolcano(
#'     object = res_small,
#'     sigPointColor = c(
#'         upregulated = "purple",
#'         downregulated = "orange"
#'     )
#' )
#'
#' # Label DEGs with a single color.
#' plotVolcano(res_small, sigPointColor = "purple")
#'
#' # Directional support (up or down).
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
#' # Label genes manually.
#' # Note that either gene IDs or names (symbols) are supported.
#' plotVolcano(
#'     object = res_small,
#'     genes = geneIDs,
#'     gene2symbol = gene2symbol
#' )
#' plotVolcano(
#'     object = res_small,
#'     genes = geneNames,
#'     gene2symbol = gene2symbol
#' )
NULL



# FIXME Improve the formals.
#' @rdname plotVolcano
#' @export
setMethod(
    "plotVolcano",
    signature("DESeqAnalysis"),
    function(object, results, ...) {
        do.call(
            what = plotVolcano,
            args = list(
                object = .matchResults(object, results),
                gene2symbol = gene2symbol(object@data),
                ...
            )
        )
    }
)



#' @rdname plotVolcano
#' @export
setMethod(
    "plotVolcano",
    signature("DESeqResults"),
    function(
        object,
        alpha = NULL,
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
        return = c("ggplot", "DataFrame")
    ) {
        validObject(object)
        if (is.null(alpha)) {
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

        # Placeholder variable for matching the LFC column.
        lfcCol <- "log2FoldChange"
        negLogTestCol <- camel(paste("neg", "log10", testCol))

        data <- object %>%
            as("tbl_df") %>%
            camel() %>%
            # Select columns used for plots.
            select(!!!syms(c("rowname", "baseMean", lfcCol, testCol))) %>%
            # Remove genes with NA adjusted P values.
            filter(!is.na(!!sym(testCol))) %>%
            # Remove genes with zero counts.
            filter(!!sym("baseMean") > 0L) %>%
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

        # Apply directional filtering, if desired.
        if (direction == "up") {
            data <- filter(data, !!sym(lfcCol) > 0L)
        } else if (direction == "down") {
            data <- filter(data, !!sym(lfcCol) < 0L)
        }

        # Early return the data, if desired.
        if (return == "DataFrame") {
            return(as(data, "DataFrame"))
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
                subtitle = paste("alpha", "<", alpha),
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
        # Get the genes to visualize when `ntop` is declared.
        if (ntop > 0L) {
            assert_is_subset(
                x = c("rowname", "rank"),
                y = colnames(data)
            )
            # Double check that data is arranged by `rank` column.
            assert_are_identical(
                x = data[["rank"]],
                y = sort(data[["rank"]])
            )
            # Since we know the data is arranged by rank, simply take the head.
            genes <- head(data[["rowname"]], n = ntop)
        }

        # Visualize specific genes on the plot, if desired.
        if (!is.null(genes)) {
            validObject(gene2symbol)
            assertFormalGene2symbol(
                object = object,
                genes = genes,
                gene2symbol = gene2symbol
            )
            # Map the user-defined `genes` to `gene2symbol` rownames.
            # We're using this to match back to the `DESeqResults` object.
            rownames <- mapGenesToRownames(
                object = gene2symbol,
                genes = genes
            )
            # Prepare the label data tibble.
            labelData <- data %>%
                .[match(x = rownames, table = .[["rowname"]]), ] %>%
                left_join(as(gene2symbol, "tbl_df"), by = "rowname")
            p <- p +
                basejump_geom_label_repel(
                    data = labelData,
                    mapping = aes(
                        x = !!sym(lfcCol),
                        y = !!sym(negLogTestCol),
                        label = !!sym("geneName")
                    )
                )
        }

        # Return ---------------------------------------------------------------
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
)
