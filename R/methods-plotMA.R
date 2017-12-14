#' Plot Mean Average
#'
#' @rdname plotMA
#' @name plotMA
#'
#' @family Differential Expression Plots
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#' @inheritParams plotHeatmap
#'
#' @importFrom BiocGenerics plotMA
#'
#' @param alpha Alpha level cutoff (Adjusted P value).
#' @param pointColor Default point color for the plot.
#' @param sigPointColor Color for points corresponding to significant genes that
#'   have passed alpha level cutoffs.
#' @param labelColor Text label color.
#' @param title *Optional*. Plot title.
#'
#' @return [ggplot].
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
#' load(system.file(
#'     file.path("extdata", "res.rda"),
#'     package = "bcbioRNASeq"))
#'
#' gene2symbol <- gene2symbol(bcb)
#'
#' ensgene <- rownames(res) %>% head(n = 4)
#' print(ensgene)
#'
#' # DESeqResults
#' plotMA(res, genes = ensgene, gene2symbol = TRUE)
#'
#' # Use a stashed gene2symbol data.frame
#' plotMA(res, genes = ensgene, gene2symbol = gene2symbol)
#'
#' # Label the Ensembl gene identifiers instead of symbols
#' plotMA(res, genes = ensgene, gene2symbol = FALSE)
#'
#' # data.frame
#' df <- as.data.frame(res)
#' plotMA(df, alpha = 0.01)
NULL



# Constructors =================================================================
#' @importFrom basejump annotable camel checkGene2symbol detectOrganism
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
    gene2symbol = TRUE,
    pointColor = "darkgray",
    sigPointColor = "purple",
    labelColor = "black",
    title = NULL) {
    data <- object %>%
        rownames_to_column("ensgene") %>%
        as_tibble() %>%
        camel(strict = FALSE) %>%
        filter(!is.na(.data[["padj"]]))
    p <- ggplot(
        data,
        mapping = aes_(
            x = ~baseMean,
            y = ~log2FoldChange,
            color = ~padj < alpha)
    ) +
        geom_point(size = 0.8) +
        scale_x_log10() +
        annotation_logticks(sides = "b") +
        guides(color = FALSE) +
        labs(title = title,
             x = "mean expression across all samples",
             y = "log2 fold change")
    if (!is.null(pointColor) & !is.null(sigPointColor)) {
        # `FALSE`: Genes that don't pass alpha
        # `TRUE`: Significant genes that do pass alpha
        p <- p +
            scale_color_manual(
                values = c("FALSE" = pointColor,
                           "TRUE" = sigPointColor))
    }
    if (!is.null(genes)) {
        if (isTRUE(gene2symbol)) {
            organism <- pull(data, "ensgene") %>%
                .[[1]] %>%
                detectOrganism()
            gene2symbol <- annotable(organism, format = "gene2symbol")
        }
        if (is.data.frame(gene2symbol)) {
            labelCol <- "symbol"
            checkGene2symbol(gene2symbol)
            data <- left_join(data, gene2symbol, by = "ensgene")
        } else {
            labelCol <- "ensgene"
        }
        labels <- data %>%
            .[.[["ensgene"]] %in% genes, , drop = FALSE]
        if (!nrow(labels)) {
            stop("Failed to label any gene identifiers")
        }
        p <- p +
            geom_text_repel(
                data = labels,
                aes_string(x = "baseMean",
                     y = "log2FoldChange",
                     label = labelCol),
                arrow = arrow(length = unit(0.01, "npc")),
                box.padding = unit(0.5, "lines"),
                color = labelColor,
                fontface = "bold",
                force = 1,
                point.padding = unit(0.75, "lines"),
                segment.color = labelColor,
                segment.size = 0.5,
                show.legend = FALSE,
                size = 4)
    }
    p
}



# Methods ======================================================================
#' @rdname plotMA
#' @importFrom ggplot2 scale_color_manual
#' @importFrom S4Vectors metadata
#' @export
setMethod(
    "plotMA",
    signature("DESeqResults"),
    function(
        object,
        genes = NULL,
        gene2symbol = TRUE,
        pointColor = "darkgray",
        sigPointColor = "red",
        labelColor = "black",
        title = TRUE) {
        results <- as.data.frame(object)
        alpha <- metadata(object)[["alpha"]]
        if (isTRUE(title)) {
            title <- .resContrastName(object)
        }
        .plotMA(
            object = results,
            alpha = alpha,
            genes = genes,
            gene2symbol = gene2symbol,
            pointColor = pointColor,
            sigPointColor = sigPointColor,
            labelColor = labelColor,
            title = title
        )
    })



#' @rdname plotMA
#' @export
setMethod(
    "plotMA",
    signature("data.frame"),
    .plotMA)
