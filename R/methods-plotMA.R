#' Plot Mean Average
#'
#' @name plotMA
#' @family Differential Expression Functions
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @importFrom BiocGenerics plotMA
#'
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @examples
#' gene2symbol <- gene2symbol(bcb_small)
#'
#' # DESeqResults ====
#' # Color DEGs in each direction separately
#' plotMA(
#'     object = res_small,
#'     sigPointColor = c(
#'         upregulated = "purple",
#'         downregulated = "orange"
#'     )
#' )
#'
#' # Label DEGs with a single color
#' plotMA(res_small, sigPointColor = "purple")
#'
#' # Directional support
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
#' # Label genes manually
#' plotMA(
#'     object = res_small,
#'     genes = head(rownames(res_small)),
#'     gene2symbol = gene2symbol
#' )
NULL



# Methods ======================================================================
#' @rdname plotMA
#' @export
setMethod(
    "plotMA",
    signature("DESeqResults"),
    function(
        object,
        genes = NULL,
        gene2symbol = NULL,
        ntop = 0L,
        direction = c("both", "up", "down"),
        pointColor = "gray50",
        sigPointColor = c(
            upregulated = "purple",
            downregulated = "orange"
        ),
        title = TRUE,
        return = c("ggplot", "data.frame")
    ) {
        validObject(object)
        assertFormalGene2symbol(object, genes, gene2symbol)
        direction <- match.arg(direction)
        assert_all_are_non_negative(ntop)
        assert_is_a_string(pointColor)
        assert_is_character(sigPointColor)
        if (is_a_string(sigPointColor)) {
            sigPointColor <- c(
                "upregulated" = sigPointColor,
                "downregulated" = sigPointColor
            )
        }
        assert_is_of_length(sigPointColor, 2L)
        assert_is_a_bool(title)
        return <- match.arg(return)

        # Title
        if (isTRUE(title)) {
            title <- contrastName(object)
        } else {
            title <- NULL
        }

        # Get alpha cutoff automatically (for coloring)
        alpha <- metadata(object)[["alpha"]]
        assert_is_a_number(alpha)
        assert_all_are_in_left_open_range(alpha, 0L, 1L)

        data <- object %>%
            as.data.frame() %>%
            rownames_to_column("geneID") %>%
            as_tibble() %>%
            camel() %>%
            mutate(rankScore = abs(!!sym("log2FoldChange"))) %>%
            arrange(desc(!!sym("rankScore"))) %>%
            mutate(rank = row_number()) %>%
            .degColors(alpha = alpha)

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

        p <- ggplot(
            data = data,
            mapping = aes_string(
                x = "baseMean",
                y = "log2FoldChange",
                color = "color"
            )
        ) +
            geom_hline(
                yintercept = 0L,
                size = 0.5,
                color = pointColor
            ) +
            geom_point(size = 1L) +
            scale_x_log10() +
            annotation_logticks(sides = "b") +
            guides(color = FALSE) +
            labs(
                title = title,
                x = "mean expression across all samples",
                y = "log2 fold change"
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
        labelData <- NULL
        if (is.null(genes) && is_positive(ntop)) {
            genes <- data[1L:ntop, "geneID", drop = TRUE]
        }
        if (is.character(genes)) {
            assert_is_subset(genes, data[["geneID"]])
            labelData <- data[data[["geneID"]] %in% genes, ]
            p <- p +
                geomLabel(
                    data = labelData,
                    mapping = aes_string(
                        x = "baseMean",
                        y = "log2FoldChange",
                        label = labelCol
                    )
                )
        }

        p
    }
)
