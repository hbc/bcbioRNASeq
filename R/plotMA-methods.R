# TODO Don't export DESeqResults method for plotMA. Consider only exporting
# DESeqAnalysis method.

# TODO Add a `results = "all"` plot mode?



#' @name plotMA
#' @inherit BiocGenerics::plotMA
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @details
#' An MA plot is an application of a Bland–Altman plot for visual representation
#' of genomic data. The plot visualizes the differences between measurements
#' taken in two samples, by transforming the data onto M (log ratio) and A
#' (mean average) scales, then plotting these values.
#'
#' @note We are not allowing post hoc `alpha` or `lfcThreshold` cutoffs here.
#'
#' @inheritParams params
#' @inheritParams basejump::params
#'
#' @return `ggplot`.
#'
#' @seealso [DESeq2::plotMA()].
#'
#' @examples
#' data(deseq)
#'
#' ## DESeqAnalysis ====
#' ## This is the current recommended method.
#' object <- deseq
#' print(object)
#'
#' dds <- as(object, "DESeqDataSet")
#' g2s <- Gene2Symbol(dds)
#' geneIDs <- head(g2s[["geneID"]])
#' print(geneIDs)
#' geneNames <- head(g2s[["geneName"]])
#' print(geneNames)
#'
#' plotMA(object)
#'
#' ## Customize the colors.
#' plotMA(
#'     object = object,
#'     pointColor = "black",
#'     sigPointColor = "purple"
#' )
#' plotMA(
#'     object = object,
#'     sigPointColor = c(
#'         upregulated = "green",
#'         downregulated = "red"
#'     )
#' )
#'
#' ## Directional support (up or down).
#' plotMA(object, direction = "up", ntop = 5L)
#' plotMA(object, direction = "down", ntop = 5L)
#'
#' ## Label genes manually.
#' ## Note that either gene IDs or names (symbols) are supported.
#' plotMA(object = object, genes = geneIDs)
#' plotMA(object = object, genes = geneNames)
NULL



#' @importFrom BiocGenerics plotMA
#' @aliases NULL
#' @export
BiocGenerics::plotMA



# DESeqResults =================================================================
plotMA.DESeqResults <-  # nolint
    function(
        object,
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
        return = c("ggplot", "DataFrame")
    ) {
        validObject(object)
        alpha <- metadata(object)[["alpha"]]
        assertIsAlpha(alpha)
        lfcThreshold <- metadata(object)[["lfcThreshold"]]
        assert_is_a_number(lfcThreshold)
        assert_all_are_non_negative(lfcThreshold)
        assert_is_any_of(genes, c("character", "NULL"))
        assert_is_any_of(gene2symbol, c("Gene2Symbol", "NULL"))
        direction <- match.arg(direction)
        assert_is_a_number(ntop)
        assert_all_are_non_negative(ntop)
        if (!is.null(genes) && ntop > 0L) {
            stop("Specify either `genes` or `ntop`.", call. = FALSE)
        }
        assert_is_a_string(pointColor)
        assert_is_character(sigPointColor)
        if (is_a_string(sigPointColor)) {
            sigPointColor <- c(
                upregulated = sigPointColor,
                downregulated = sigPointColor
            )
        }
        assert_is_of_length(sigPointColor, n = 2L)
        return <- match.arg(return)

        # Check to see if we should use `sval` column instead of `padj`.
        if ("svalue" %in% names(object)) {
            testCol <- "svalue"  # nocov
        } else {
            testCol <- "padj"
        }

        # Placeholder variable for matching the LFC column.
        lfcCol <- "log2FoldChange"

        data <- object %>%
            as_tibble(rownames = "rowname") %>%
            camel() %>%
            # Remove genes with very low expression.
            filter(!!sym("baseMean") >= 1L) %>%
            mutate(rankScore = abs(!!sym("log2FoldChange"))) %>%
            arrange(desc(!!sym("rankScore"))) %>%
            mutate(rank = row_number()) %>%
            .addIsDECol(
                testCol = testCol,
                alpha = alpha,
                lfcCol = lfcCol,
                lfcThreshold = lfcThreshold
            )
        assert_is_subset(
            x = c(
                "rowname",
                "baseMean",
                lfcCol,
                testCol,
                "rankScore",
                "rank",
                "isDE"
            ),
            y = colnames(data)
        )

        # Apply directional filtering, if desired.
        if (direction == "up") {
            data <- filter(data, !!sym(lfcCol) > 0L)
        } else if (direction == "down") {
            data <- filter(data, !!sym(lfcCol) < 0L)
        }

        # Check for no genes passing cutoffs and early return.
        if (!nrow(data)) {
            warning("No genes passed cutoffs.")
            return(invisible())
        }

        # Early return the data, if desired.
        if (return == "DataFrame") {
            return(as(data, "DataFrame"))
        }

        # MA plot --------------------------------------------------------------
        log10BaseMean <- log10(data[["baseMean"]])
        floor <- min(floor(log10BaseMean))
        ceiling <- max(ceiling(log10BaseMean))
        xBreaks <- 10L ^ seq(from = floor, to = ceiling, by = 1L)

        p <- ggplot(
            data = data,
            mapping = aes(
                x = !!sym("baseMean"),
                y = !!sym(lfcCol),
                color = !!sym("isDE")
            )
        ) +
            geom_hline(
                yintercept = 0L,
                size = 0.5,
                color = pointColor
            ) +
            geom_point(size = 1L) +
            scale_x_continuous(
                breaks = xBreaks,
                limits = c(1L, NA),
                trans = "log10"
            ) +
            scale_y_continuous(breaks = pretty_breaks()) +
            annotation_logticks(sides = "b") +
            guides(color = FALSE) +
            labs(
                title = contrastName(object),
                subtitle = paste("alpha", "<", alpha),
                x = "mean expression across all samples",
                y = "log2 fold change"
            )

        # Color the significant points.
        # Note that we're using direction-specific coloring by default.
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
                        x = !!sym("baseMean"),
                        y = !!sym(lfcCol),
                        label = !!sym("geneName")
                    )
                )
        }

        # Return ---------------------------------------------------------------
        p
    }



#' @rdname plotMA
#' @export
setMethod(
    f = "plotMA",
    signature = signature("DESeqResults"),
    definition = plotMA.DESeqResults
)



# DESeqAnalysis ================================================================
plotMA.DESeqAnalysis <-  # nolint
    function(
        object,
        results = 1L,
        lfcShrink = TRUE
    ) {
        validObject(object)
        do.call(
            what = plotMA,
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
f1 <- formals(plotMA.DESeqAnalysis)
f2 <- formals(plotMA.DESeqResults)
f2 <- f2[setdiff(names(f2), c(names(f1), "gene2symbol"))]
f <- c(f1, f2)
formals(plotMA.DESeqAnalysis) <- f



#' @rdname plotMA
#' @export
setMethod(
    f = "plotMA",
    signature = signature("DESeqAnalysis"),
    definition = plotMA.DESeqAnalysis
)



# Aliases ======================================================================
#' @rdname plotMA
#' @export
plotMeanAverage <- function(...) {
    # This function is soft deprecated.
    plotMA(...)
}
