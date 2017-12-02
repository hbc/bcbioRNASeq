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
#' @param countsAxisLabel Text label of counts axis.
#' @param returnList Return the plotlist used to generate the paneled,
#'   multi-gene plot with [cowplot::plot_grid()].
#' @param metadata Sample metadata [data.frame].
#'
#' @return
#' - `returnList = FALSE`: [cowplot::plot_grid()] graphical output.
#' - `returnList = TRUE`: [list] of per gene [ggplot] objects.
#'
#' @seealso [DESeq2::plotCounts()].
#'
#' @examples
#' bcb <- examples[["bcb"]]
#'
#' # Gene symbols
#' symbol <- rowData(bcb)[["symbol"]][1:3]
#' print(symbol)
#' plotGene(bcb, gene = symbol, format = "symbol")
#' plotGene(
#'     bcb,
#'     gene = symbol,
#'     format = "symbol",
#'     interestingGroups = "sampleName",
#'     color = NULL)
#'
#' # Gene identifiers
#' \dontrun{
#' ensgene <- rowData(bcb)[["ensgene"]][1:4]
#' print(ensgene)
#' plotGene(bcb, gene = ensgene, format = "ensgene")
#' }
NULL



# Constructors ====
#' Plot Gene Constructor
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @importFrom basejump uniteInterestingGroups
#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 aes_string element_text expand_limits geom_point ggplot
#'   labs theme
#' @importFrom viridis scale_color_viridis
#'
#' @param object Counts matrix.
#' @param gene Gene identifiers, as a named character vector. `ensgene` is
#'   provided as the value and `symbol` as the name. This gets defined in the S4
#'   method (see below).
#' @param metadata Sample metadata [data.frame].
#'
#' @return [ggplot].
.plotGene <- function(
    object,
    gene,
    metadata,
    interestingGroups = "sampleName",
    color = viridis::scale_color_viridis(discrete = TRUE),
    countsAxisLabel = "counts",
    returnList = FALSE) {
    metadata <- metadata %>%
        as.data.frame() %>%
        uniteInterestingGroups(interestingGroups)
    plots <- lapply(seq_along(gene), function(a) {
        ensgene <- gene[[a]]
        symbol <- names(gene)[[a]]
        df <- data.frame(
            x = colnames(object),
            y = object[ensgene, ],
            interestingGroups = metadata[["interestingGroups"]])
        p <- ggplot(
            df,
            mapping = aes_string(
                x = "x",
                y = "y",
                color = "interestingGroups")
        ) +
            geom_point(size = 4) +
            theme(
                axis.text.x = element_text(angle = 90)) +
            labs(title = symbol,
                 x = "sample",
                 y = countsAxisLabel,
                 color = paste(interestingGroups, collapse = ":\n")) +
            expand_limits(y = 0)
        if (!is.null(color)) {
            p <- p + color
        }
        p
    })
    if (isTRUE(returnList)) {
        plots
    } else {
        plot_grid(plotlist = plots, labels = "AUTO")
    }
}



# Methods ====
#' @rdname plotGene
#' @importFrom viridis scale_color_viridis
#' @export
setMethod(
    "plotGene",
    signature("bcbioRNASeq"),
    function(
        object,
        gene,
        interestingGroups,
        normalized = "tpm",
        format = "ensgene",
        color = viridis::scale_color_viridis(discrete = TRUE),
        returnList = FALSE) {
        if (!format %in% c("ensgene", "symbol")) {
            stop("Unsupported gene identifier format", call. = FALSE)
        }
        supportedAssay <- c("tpm", "tmm", "rlog", "vst")
        if (!normalized %in% supportedAssay) {
            stop(paste(
                "'normalized' argument requires:",
                toString(supportedAssay)
            ), call. = FALSE)
        }
        countsAxisLabel <- normalized
        if (missing(interestingGroups)) {
            interestingGroups <- basejump::interestingGroups(object)
        }

        counts <- counts(object, normalized = normalized)
        metadata <- sampleMetadata(object)

        # Match unique gene identifier with name (gene symbol) using the
        # internally stored Ensembl annotations saved in the run object
        gene2symbol <- gene2symbol(object)

        # Detect missing genes. This also handles `format` mismatch.
        if (!all(gene %in% gene2symbol[[format]])) {
            stop(paste(
                "Missing genes:",
                toString(setdiff(gene, gene2symbol[[format]]))
            ), call. = FALSE)
        }

        # Now safe to prepare the named gene character vector. Here we're
        # passing in the `ensgene` as the value and `symbol` as the name. This
        # works with the constructor function to match the counts matrix by
        # the ensgene, then use the symbol as the name.
        match <- match(x = gene, table = gene2symbol[[format]])
        gene2symbol <- gene2symbol[match, , drop = FALSE]
        ensgene <- gene2symbol[["ensgene"]]
        names(ensgene) <- gene2symbol[["symbol"]]

        .plotGene(
            object = counts,
            gene = ensgene,
            metadata = metadata,
            interestingGroups = interestingGroups,
            color = color,
            countsAxisLabel = countsAxisLabel,
            returnList = returnList)
    })



#' @rdname plotGene
#' @export
setMethod(
    "plotGene",
    signature("matrix"),
    .plotGene)
