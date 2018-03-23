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
#' @param sigPointColor Color for points corresponding to significant genes that
#'   have passed alpha level cutoffs.
#' @param labelColor Text label color.
#'
#' @return `ggplot`.
#'
#' @examples
#' load(system.file("extdata/bcb_small.rda", package = "bcbioRNASeq"))
#' load(system.file("extdata/res_small.rda", package = "bcbioRNASeq"))
#'
#' # DESeqResults ====
#' plotMA(res_small)
#' plotMA(
#'     object = res_small,
#'     genes = head(rownames(res_small), 4L),
#'     gene2symbol = gene2symbol(bcb_small)
#' )
NULL



# Constructors =================================================================
.plotMA <- function(
    object,
    genes = NULL,
    gene2symbol = NULL,
    pointColor = "darkgray",
    sigPointColor = "purple",
    labelColor = "black",
    title = TRUE
) {
    validObject(object)
    assertFormalGene2symbol(object, genes, gene2symbol)
    assert_is_a_string(pointColor)
    assert_is_a_string(sigPointColor)
    assert_is_a_string(labelColor)
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
        .[!is.na(.[["padj"]]), , drop = FALSE]

    p <- ggplot(
        data = data,
        mapping = aes_(
            x = ~baseMean,
            y = ~log2FoldChange,
            color = ~padj < alpha
        )
    ) +
        geom_hline(
            yintercept = 0L,
            size = 1L,
            color = sigPointColor,
            alpha = 0.25
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
    .plotMA
)
