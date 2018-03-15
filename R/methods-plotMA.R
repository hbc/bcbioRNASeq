#' Plot Mean Average
#'
#' @name plotMA
#' @family Differential Expression Functions
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @inheritParams general
#' @inheritParams plotHeatmap
#'
#' @importFrom BiocGenerics plotMA
#'
#' @param pointColor Default point color for the plot.
#' @param sigPointColor Color for points corresponding to significant genes that
#'   have passed alpha level cutoffs.
#' @param labelColor Text label color.
#'
#' @return [ggplot].
#'
#' @examples
#' # DESeqResults ====
#' plotMA(res_small)
#' plotMA(
#'     object = res_small,
#'     genes = head(rownames(res_small), 4L),
#'     gene2symbol = gene2symbol(bcb_small)
#' )
NULL



# Constructors =================================================================
#' @importFrom basejump camel detectOrganism
#' @importFrom dplyr filter pull
#' @importFrom ggplot2 aes_ annotation_logticks geom_point ggtitle guides labs
#'   scale_color_manual scale_x_log10
#' @importFrom ggrepel geom_text_repel
#' @importFrom grid arrow unit
#' @importFrom tibble as_tibble rownames_to_column
.plotMA <- function(
    object,
    alpha = 0.01,
    genes = NULL,
    gene2symbol = NULL,
    pointColor = "darkgray",
    sigPointColor = "purple",
    labelColor = "black",
    title = NULL
) {
    assert_is_data.frame(object)
    assert_is_a_number(alpha)
    assert_all_are_positive(alpha)
    assertFormalGene2symbol(object, genes, gene2symbol)
    assert_is_a_string(pointColor)
    assert_is_a_string(sigPointColor)
    assert_is_a_string(labelColor)
    assertIsAStringOrNULL(title)

    data <- object %>%
        rownames_to_column("geneID") %>%
        as_tibble() %>%
        camel(strict = FALSE) %>%
        filter(!is.na(.data[["padj"]]))

    p <- ggplot(
        data = data,
        mapping = aes_(
            x = ~baseMean,
            y = ~log2FoldChange,
            color = ~padj < alpha)
    ) +
        geom_point(size = 0.8) +
        scale_x_log10() +
        annotation_logticks(sides = "b") +
        guides(color = FALSE) +
        labs(
            title = title,
            x = "mean expression across all samples",
            y = "log2 fold change"
        )

    if (!is.null(pointColor) && !is.null(sigPointColor)) {
        # `FALSE`: Genes that don't pass alpha
        # `TRUE`: Significant genes that do pass alpha
        p <- p +
            scale_color_manual(
                values = c("FALSE" = pointColor, "TRUE" = sigPointColor)
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
            geom_text_repel(
                data = labels,
                aes_string(
                    x = "baseMean",
                    y = "log2FoldChange",
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

    p
}



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
        pointColor = "darkgray",
        sigPointColor = "red",
        labelColor = "black",
        title = TRUE
    ) {
        # Passthrough: genes, gene2symbol, pointColor, sigPointColor, labelColor
        if (isTRUE(title)) {
            title <- contrastName(object)
        } else if (!is_a_string(title)) {
            title <- NULL
        }
        .plotMA(
            object = as.data.frame(object),
            alpha = metadata(object)[["alpha"]],
            genes = genes,
            gene2symbol = gene2symbol,
            pointColor = pointColor,
            sigPointColor = sigPointColor,
            labelColor = labelColor,
            title = title
        )
    }
)
