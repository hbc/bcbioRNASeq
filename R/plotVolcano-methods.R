#' @name plotVolcano
#' @inherit basejump::plotVolcano
#' @author Michael Steinbaugh, John Hutchinson, Lorena Pantano
#'
#' @inheritParams params
#' @inheritParams basejump::params
#' @param ylim `scalar numeric`. Upper boundary limit for y-axis. Helps preserve
#'   dynamic range for gene sets containing highly significant P values (e.g.
#'   `1e-100`).
#' @param histograms `boolean`. Show LFC and P value histograms.
#'
#' @seealso This method is an updated variant of
#'   `CHBUtils::volcano_density_plot()`.
#'
#' @examples
#' data(deseq)
#'
#' object <- deseq
#' print(object)
#'
#' ## Get genes from DESeqDataSet.
#' dds <- as(object, "DESeqDataSet")
#' g2s <- Gene2Symbol(dds)
#' geneIDs <- head(g2s[["geneID"]])
#' print(geneIDs)
#' geneNames <- head(g2s[["geneName"]])
#' print(geneNames)
#'
#' plotVolcano(object)
#'
#' ## Customize the colors.
#' plotVolcano(
#'     object = object,
#'     pointColor = "black",
#'     sigPointColor = "purple"
#' )
#' plotVolcano(
#'     object = object,
#'     sigPointColor = c(
#'         upregulated = "green",
#'         downregulated = "red"
#'     )
#' )
#'
#' ## Directional support (up or down).
#' plotVolcano(
#'     object = object,
#'     direction = "up",
#'     ntop = 5L
#' )
#' plotVolcano(
#'     object = object,
#'     direction = "down",
#'     ntop = 5L
#' )
#'
#' ## Label genes manually.
#' ## Note that either gene IDs or names (symbols) are supported.
#' plotVolcano(object, genes = geneIDs)
#' plotVolcano(object, genes = geneNames)
NULL



#' @importFrom basejump plotVolcano
#' @aliases NULL
#' @export
basejump::plotVolcano



# DESeqResults =================================================================
plotVolcano.DESeqResults <-  # nolint
    function(
        object,
        ylim = 1e-10,
        genes = NULL,
        gene2symbol = NULL,
        ntop = 0L,
        direction = c("both", "up", "down"),
        pointColor = getOption("basejump.color", "gray50"),
        sigPointColor = getOption(
            "basejump.point.color",
            c(
                upregulated = "purple",
                downregulated = "orange"
            )
        ),
        histograms = FALSE,
        return = c("ggplot", "DataFrame")
    ) {
        validObject(object)
        alpha <- metadata(object)[["alpha"]]
        assertIsAlpha(alpha)
        lfcThreshold <- metadata(object)[["lfcThreshold"]]
        assert_is_a_number(lfcThreshold)
        assert_is_a_number(ylim)
        assert_all_are_in_range(
            x = ylim,
            lower = 1e-100,
            upper = 1e-3
        )
        assertIsImplicitInteger(ntop)
        assert_all_are_non_negative(c(lfcThreshold, ntop))
        direction <- match.arg(direction)
        assert_is_a_string(pointColor)
        assert_is_character(sigPointColor)
        if (is_a_string(sigPointColor)) {
            sigPointColor <- c(
                upregulated = sigPointColor,
                downregulated = sigPointColor
            )
        }
        assert_is_of_length(sigPointColor, n = 2L)
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
            as_tibble(rownames = "rowname") %>%
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
            assertFormalGene2Symbol(
                x = object,
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
                    x = 0L, y = 0.2,
                    width = 1L, height = 0.8
                ) +
                draw_plot(
                    plot = lfcHist,
                    x = 0L, y = 0L,
                    width = 0.5, height = 0.2
                ) +
                draw_plot(
                    plot = pvalueHist,
                    x = 0.5, y = 0L,
                    width = 0.5, height = 0.2
                )
        } else {
            p
        }
    }



#' @rdname plotVolcano
#' @export
setMethod(
    "plotVolcano",
    signature("DESeqResults"),
    definition = plotVolcano.DESeqResults
)



# DESeqAnalysis ================================================================
plotVolcano.DESeqAnalysis <-  # nolint
    function(
        object,
        results = 1L,
        lfcShrink = TRUE
    ) {
        validObject(object)
        do.call(
            what = plotVolcano,
            args = matchArgsToDoCall(
                args = list(
                    object = .matchResults(
                        object = object,
                        results = results,
                        lfcShrink = lfcShrink
                    ),
                    genes = genes,
                    gene2symbol = Gene2Symbol(slot(object, "data"))
                ),
                removeFormals = c("results", "lfcShrink")
            )
        )
    }
f1 <- formals(plotVolcano.DESeqAnalysis)
f2 <- formals(plotVolcano.DESeqResults)
f2 <- f2[setdiff(names(f2), c(names(f1), "gene2symbol"))]
f <- c(f1, f2)
formals(plotVolcano.DESeqAnalysis) <- f



#' @rdname plotVolcano
#' @export
setMethod(
    f = "plotVolcano",
    signature = signature("DESeqAnalysis"),
    definition = plotVolcano.DESeqAnalysis
)
