#' Plot Individual Genes
#'
#' @rdname plotGene
#' @name plotGene
#' @family Quality Control Plots
#' @author Michael Steinbaugh
#'
#' @importFrom basejump plotGene
#'
#' @inherit plotTotalReads
#'
#' @param gene Gene identifier. Can input multiple genes as a character vector.
#' @param format Ensembl identifier format. Supports `symbol` (**default**) or
#'   `ensgene`.
#' @param normalized Normalization method. Supports `tpm` (**default**), `tmm`,
#'   `rlog`, or `vst`.
#' @param color Desired ggplot color scale. Defaults to
#'   [viridis::scale_color_viridis()]. Must supply discrete values. When set to
#'   `NULL`, the default ggplot2 color palette will be used. If manual color
#'   definitions are desired, we recommend using
#'   [ggplot2::scale_color_manual()].
#'
#' @return
#' - `return = FALSE`: [cowplot::plot_grid()] graphical output.
#' - `return = TRUE`: [list] of per gene [ggplot] objects.
#'
#' @seealso [DESeq2::plotCounts()].
#'
#' @examples
#' # Gene symbols
#' symbol <- rowData(bcb)[["symbol"]][1:3]
#' print(symbol)
#' plotGene(bcb, gene = symbol, format = "symbol")
#'
#' # Gene identifiers
#' \dontrun{
#' ensgene <- rownames(bcb)[1:3]
#' print(ensgene)
#' plotGene(bcb, gene = ensgene, format = "ensgene")
#' }
NULL



# Constructors ====
#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 aes_string element_text expand_limits geom_point ggplot
#'   labs theme
#' @importFrom viridis scale_color_viridis
.plotGene <- function(
    object,
    gene,
    metadata,
    interestingGroups = "sampleName",
    color = scale_color_viridis(discrete = TRUE),
    return = FALSE) {
    interestingGroups <- .checkInterestingGroups(object, interestingGroups)
    metadata <- as.data.frame(metadata)
    plots <- lapply(seq_along(gene), function(a) {
        ensgene <- gene[[a]]
        symbol <- names(gene)[[a]]
        df <- data.frame(
            x = colnames(counts),
            y = counts[ensgene, ],
            color = metadata[[interestingGroups]])
        p <- ggplot(
            df,
            mapping = aes_string(
                x = "x",
                y = "y",
                color = "color")
        ) +
            geom_point(size = 4) +
            theme(
                axis.text.x = element_text(angle = 90)) +
            labs(title = symbol,
                 x = "sample",
                 y = "counts",
                 color = interestingGroups) +
            expand_limits(y = 0)
        if (!is.null(color)) {
            p <- p + color
        }
        p
    })
    if (isTRUE(return)) {
        plots
    } else {
        plot_grid(plotlist = plots, labels = "AUTO")
    }
}



# Methods ====
#' @rdname plotGene
#' @importFrom S4Vectors metadata
#' @importFrom viridis scale_color_viridis
#' @export
setMethod(
    "plotGene",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        normalized = "tpm",
        gene,
        format = "symbol",
        color = scale_color_viridis(discrete = TRUE)) {
        if (!format %in% c("ensgene", "symbol")) {
            stop("Unsupported gene identifier format", call. = FALSE)
        }
        if (is.logical(normalized)) {
            warning(paste(
                "Explicit format of normalized counts format is recommended"
            ), call. = FALSE)
        }
        if (missing(interestingGroups)) {
            interestingGroups <-
                metadata(object)[["interestingGroups"]][[1]]
        }

        counts <- counts(object, normalized = normalized)
        metadata <- colData(object)

        # Match unique gene identifier with name (gene symbol) using the
        # internally stored Ensembl annotations saved in the run object
        gene2symbol <- metadata(object) %>%
            .[["annotable"]] %>%
            .[.[[format]] %in% gene, , drop = FALSE] %>%
            .[, c("symbol", "ensgene")]
        gene <- gene2symbol[["ensgene"]]
        names(gene) <- gene2symbol[["symbol"]]

        .plotGene(
            object = counts,
            gene = gene,
            metadata = metadata,
            interestingGroups = interestingGroups,
            color = color)
    })
