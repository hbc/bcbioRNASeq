# FIXME Include the alpha information in the plot.
# FIXME Use `mapGenesToRownames()` to handle `genes` input.
# FIXME Improve gene2symbol mapping support (see plotVolcano.)
# FIXME `useMapGenesToRownames()` for `DataFrame` method.



#' MA Plot
#'
#' An MA plot is an application of a Bland–Altman plot for visual representation
#' of genomic data. The plot visualizes the differences between measurements
#' taken in two samples, by transforming the data onto M (log ratio) and A
#' (mean average) scales, then plotting these values.
#'
#' @name plotMA
#' @family Differential Expression Functions
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @importFrom BiocGenerics plotMA
#' @export
#'
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @seealso [DESeq2::plotMA()].
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
#' plotMA(
#'     object = res_small,
#'     sigPointColor = c(
#'         upregulated = "purple",
#'         downregulated = "orange"
#'     )
#' )
#'
#' # Label DEGs with a single color.
#' plotMA(res_small, sigPointColor = "purple")
#'
#' # Directional support (up or down).
#' plotMA(
#'     object = res_small,
#'     direction = "up",
#'     ntop = 5L,
#'     gene2symbol = gene2symbol
#' )
#' plotMA(
#'     object = res_small,
#'     direction = "down",
#'     ntop = 5L,
#'     gene2symbol = gene2symbol
#' )
#'
#' # Label genes manually.
#' # Note that either gene IDs or names (symbols) are supported.
#' plotMA(
#'     object = res_small,
#'     genes = geneIDs,
#'     gene2symbol = gene2symbol
#' )
#' plotMA(
#'     object = res_small,
#'     genes = geneNames,
#'     gene2symbol = gene2symbol
#' )
NULL



#' @rdname plotMA
#' @export
setMethod(
    "plotMA",
    signature("DESeqResults"),
    function(
        object,
        alpha = NULL,
        lfcThreshold = 0L,
        genes = NULL,
        gene2symbol = NULL,
        ntop = 0L,
        direction = c("both", "up", "down"),
        pointColor = "gray50",
        sigPointColor = c(
            upregulated = "purple",
            downregulated = "orange"
        ),
        return = c("ggplot", "DataFrame")
    ) {
        validObject(object)
        if (is.null(alpha)) {
            alpha <- metadata(object)[["alpha"]]
        }
        assert_is_a_number(alpha)
        assert_all_are_in_left_open_range(alpha, lower = 0L, upper = 1L)
        assert_is_a_number(lfcThreshold)
        assert_all_are_non_negative(lfcThreshold)
        assert_is_any_of(genes, c("character", "NULL"))
        assert_is_any_of(gene2symbol, c("gene2symbol", "NULL"))
        direction <- match.arg(direction)
        assert_all_are_non_negative(ntop)
        if (!is.null(genes) && !is.null(ntop)) {
            stop("Specify either `genes` or `ntop`, but not both")
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
            as("tbl_df") %>%
            camel() %>%
            # Remove genes with very low expression
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
            data <- filter(data, sym(lfcCol) > 0L)
        } else if (direction == "down") {
            data <- filter(data, sym(lfcCol) < 0L)
        }

        # Check for no genes passing cutoffs and early return.
        if (!nrow(data)) {
            warning("No genes passed cutoffs")
            return(invisible())
        }

        # Early return data frame, if desired.
        if (return == "DataFrame") {
            return(as(data, "DataFrame"))
        }

        # ggplot ---------------------------------------------------------------
        xFloor <- data[["baseMean"]] %>%
            min() %>%
            log10() %>%
            floor()
        xCeiling <- data[["baseMean"]] %>%
            max() %>%
            log10() %>%
            ceiling()
        xBreaks <- 10L ^ seq(from = xFloor, to = xCeiling, by = 1L)

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
        if (length(ntop)) {
            assert_is_subset(
                x = c("rowname", "rank"),
                y = colnames(data)
            )
            # Double check that data is arranged by `rank` column.
            assert_are_identical(
                x = seq_len(nrow(data)),
                y = data[["rank"]]
            )
            # Since we know the data is arranged by rank, simply take the head.
            genes <- head(data[["rowname"]], n = ntop)
        }

        # Visualize specific genes on the plot, if desired.
        if (length(genes)) {
            validObject(gene2symbol)
            assertFormalGene2symbol(
                object = object,
                genes = genes,
                gene2symbol = gene2symbol
            )
            # Map the user-defined `genes` to `gene2symbol` rownames.
            # We're using this to match back to the `DESeqResults` object.
            rownames <- mapGenesToRownames(object = gene2symbol, genes = genes)
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
)
