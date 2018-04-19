#' Plot Mean Average
#'
#' @name plotMA
#' @family Differential Expression Functions
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @importFrom BiocGenerics plotMA
#'
#' @inheritParams general
#' @param pointColor Default point color for the plot.
#' @param sigPointColor `character` vector containing color names for labeling
#'   upregulated and downregulated genes. Also supports a character string for
#'   labeling DEGs with the same color, regardless of direction.
#' @param labelColor Text label color.
#'
#' @return `ggplot`.
#'
#' @examples
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
#' # Label DEGs in both directions with a single color
#' plotMA(
#'     object = res_small,
#'     sigPointColor = "purple"
#' )
#'
#' # Label genes
#' genes <- head(rownames(res_small), 4L)
#' gene2symbol <- gene2symbol(bcb_small)
#' plotMA(
#'     object = res_small,
#'     genes = genes,
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
        pointColor = "gray50",
        sigPointColor = c(
            upregulated = "purple",
            downregulated = "orange"
        ),
        title = TRUE
    ) {
        validObject(object)
        assertFormalGene2symbol(object, genes, gene2symbol)
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

        # Get alpha cutoff automatically
        alpha <- metadata(object)[["alpha"]]
        assert_is_a_number(alpha)
        assert_all_are_in_left_open_range(alpha, 0L, 1L)

        # Title
        if (isTRUE(title)) {
            title <- contrastName(object)
        } else {
            title <- NULL
        }

        data <- object %>%
            as.data.frame() %>%
            rownames_to_column("geneID") %>%
            as_tibble() %>%
            camel(strict = FALSE) %>%
            .[!is.na(.[["padj"]]), , drop = FALSE] %>%
            .degColors(alpha = alpha)

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

        if (is.character(genes)) {
            if (is.data.frame(gene2symbol)) {
                labelCol <- "geneName"
                assertIsGene2symbol(gene2symbol)
                data <- left_join(data, gene2symbol, by = "geneID")
            } else {
                labelCol <- "geneID"
            }
            labels <- data %>%
                .[.[["geneID"]] %in% genes, , drop = FALSE]
            assert_is_non_empty(labels)
            p <- p +
                geomLabel(
                    data = labels,
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
