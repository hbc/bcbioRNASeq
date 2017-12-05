#' Plot Mean Average
#'
#' @rdname plotMA
#' @name plotMA
#'
#' @family Differential Expression Plots
#' @author Rory Kirchner, Michael Steinbaugh
#' @inheritParams AllGenerics
#'
#' @importFrom BiocGenerics plotMA
#'
#' @param alpha Alpha level cutoff (Adjusted P value).
#' @param genes *Optional*. Genes to label on the plot.
#' @param format Gene identifier format (`ensgene` or `symbol`) used for
#'   matching.
#' @param gene2symbol Gene to symbol [data.frame] that defines
#'   identifier mappings. Required if `genes` is defined.
#' @param color Default point color for the plot.
#' @param sigColor Color for points corresponding to significant genes that have
#'   passed alpha level cutoffs.
#' @param labelColor Text label color.
#' @param title *Optional*. Plot title.
#'
#' @return [ggplot].
#'
#' @examples
#' bcb <- examples[["bcb"]]
#' res <- examples[["res"]]
#' gene2symbol <- gene2symbol(bcb)
#'
#' # DESeqResults
#' plotMA(
#'     res,
#'     genes = "ENSMUSG00000016918",
#'     format = "ensgene",
#'     gene2symbol = gene2symbol)
#' plotMA(
#'     res,
#'     genes = "Sulf1",
#'     format = "symbol",
#'     gene2symbol = gene2symbol)
#'
#' # data.frame
#' \dontrun{
#' df <- as.data.frame(res)
#' plotMA(df, alpha = 0.01)
#' }
NULL



# Constructors ====
#' @importFrom basejump camel
#' @importFrom dplyr filter
#' @importFrom ggplot2 aes_ annotation_logticks geom_point ggtitle guides labs
#'   scale_color_manual scale_x_log10
#' @importFrom ggrepel geom_text_repel
#' @importFrom grid arrow unit
#' @importFrom tibble as_tibble rownames_to_column
.plotMA <- function(
    object,
    alpha = 0.01,
    genes = NULL,
    format = "ensgene",
    gene2symbol = NULL,
    color = "darkgray",
    sigColor = "red",
    labelColor = "black",
    title = NULL) {
    validFormat <- c("ensgene", "symbol")
    if (!format %in% validFormat) {
        stop(paste("format:", toString(validFormat)), call. = FALSE)
    }
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
    if (!is.null(color) & !is.null(sigColor)) {
        # `FALSE`: Genes that don't pass alpha
        # `TRUE`: Significant genes that do pass alpha
        p <- p +
            scale_color_manual(
                values = c("FALSE" = color,
                           "TRUE" = sigColor))
    }
    if (!is.null(genes)) {
        if (!is.data.frame(gene2symbol)) {
            stop("'gene2symbol' data.frame required when 'genes' is defined",
                 call. = FALSE)
        }
        if (!all(colnames(gene2symbol) %in% validFormat)) {
            stop(paste(
                "'gene2symbol' must contain:", toString(validFormat)
            ), call. = FALSE)
        }
        data <- left_join(data, gene2symbol, by = "ensgene")
        labels <- data %>%
            .[.[[format]] %in% genes, , drop = FALSE]
        if (!nrow(labels)) {
            stop("Failed to label any gene identifiers")
        }
        p <- p +
            geom_text_repel(
                data = labels,
                aes_string(x = "baseMean",
                     y = "log2FoldChange",
                     label = "symbol"),
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



# Methods ====
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
        format = "ensgene",
        gene2symbol = NULL,
        color = "darkgray",
        sigColor = "red",
        labelColor = "black") {
        .plotMA(
            object = as.data.frame(object),
            alpha = metadata(object)[["alpha"]],
            genes = genes,
            format = format,
            gene2symbol = gene2symbol,
            color = color,
            sigColor = sigColor,
            labelColor = labelColor,
            title = .resContrastName(object)
        )
    })



#' @rdname plotMA
#' @export
setMethod(
    "plotMA",
    signature("data.frame"),
    .plotMA)
