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
#' @param genes Gene identifier(s). Can input multiple genes as a character
#'   vector.
#' @param format Ensembl identifier format. Supports `ensgene` (**recommended**)
#'   or `symbol`.
#' @param metadata Sample metadata [data.frame].
#' @param normalized Normalization method. Supports `tpm` (**default**), `tmm`,
#'   `rlog`, or `vst`.
#' @param color Desired ggplot color scale. Defaults to
#'   [viridis::scale_color_viridis()]. Must supply discrete values. When set to
#'   `NULL`, the default ggplot2 color palette will be used. If manual color
#'   definitions are desired, we recommend using
#'   [ggplot2::scale_color_manual()].
#' @param countsAxisLabel Text label of counts axis.
#' @param headerLevel R Markdown header level. Only applies when
#'   `return = "markdown"`.
#' @param return Desired return type: `grid`, `list`, `markdown`.
#'
#' @return
#' - `grid`: [cowplot::plot_grid()] graphical output.
#' - `list`: Plot list, containing per gene [ggplot] objects.
#' - `markdown`: Tabset R Markdown output, tabbed per gene.
#'
#' @seealso [DESeq2::plotCounts()].
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
#'
#' # Gene identifiers
#' ensgene <- rownames(bcb)[1:4]
#' print(ensgene)
#' plotGene(
#'     bcb,
#'     genes = ensgene,
#'     format = "ensgene")
#'
#' # Gene symbols
#' symbol <- rowData(bcb)[["symbol"]][1:4]
#' print(symbol)
#' plotGene(
#'     bcb,
#'     genes = symbol,
#'     format = "symbol")
#'
#' # Default ggplot2 color palette
#' plotGene(
#'     bcb,
#'     genes = ensgene,
#'     format = "ensgene",
#'     interestingGroups = "sampleName",
#'     color = NULL)
NULL



# Constructors =================================================================
#' Plot Gene Constructor
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @importFrom basejump uniteInterestingGroups
#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 aes_string element_text expand_limits geom_point ggplot
#'   guides labs theme
#' @importFrom pbapply pblapply
#' @importFrom tibble tibble
#' @importFrom viridis scale_color_viridis
#'
#' @param genes Gene identifiers, as a named character vector. `ensgene` is
#'   provided as the value and `symbol` as the name. This gets defined in the S4
#'   method (see below).
#'
#' @return [ggplot].
.plotGene <- function(
    object,
    genes,
    metadata,
    interestingGroups = "sampleName",
    color = viridis::scale_color_viridis(discrete = TRUE),
    countsAxisLabel = "counts",
    headerLevel = 2,
    return = "grid") {
    metadata <- metadata %>%
        as.data.frame() %>%
        uniteInterestingGroups(interestingGroups)
    plots <- lapply(seq_along(genes), function(a) {
        ensgene <- genes[[a]]
        symbol <- names(genes)[[a]]
        data <- tibble(
            x = colnames(object),
            y = object[ensgene, ],
            interestingGroups = metadata[["interestingGroups"]])
        p <- ggplot(
            data,
            mapping = aes_string(
                x = "x",
                y = "y",
                color = "interestingGroups")
        ) +
            geom_point(size = 4) +
            theme(axis.text.x = element_text(angle = 90)) +
            labs(title = symbol,
                 x = "sample",
                 y = countsAxisLabel,
                 color = paste(interestingGroups, collapse = ":\n")) +
            expand_limits(y = 0)
        if (is(color, "ScaleDiscrete")) {
            p <- p + color
        }
        if (interestingGroups == "sampleName") {
            p <- p + guides(color = FALSE)
        }
        p
    })
    names(plots) <- genes

    # Return
    validReturn <- c("grid", "list", "markdown")
    if (return == "grid") {
        plot_grid(plotlist = plots, labels = "AUTO")
    } else if (return == "list") {
        plots
    } else if (return == "markdown") {
        pblapply(seq_along(plots), function(a) {
            mdHeader(names(genes)[[a]], level = headerLevel)
            plots[[a]]
        }) %>%
            invisible()
    } else {
        stop(paste("Valid 'return':", toString(validReturn)))
    }
}



# Methods ======================================================================
#' @rdname plotGene
#' @importFrom viridis scale_color_viridis
#' @export
setMethod(
    "plotGene",
    signature("bcbioRNASeq"),
    function(
        object,
        genes,
        format = "ensgene",
        normalized = "tpm",
        interestingGroups,
        color = viridis::scale_color_viridis(discrete = TRUE),
        headerLevel = 2,
        return = "grid") {
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
        if (!all(genes %in% gene2symbol[[format]])) {
            stop(paste(
                "Missing genes:",
                toString(setdiff(genes, gene2symbol[[format]]))
            ), call. = FALSE)
        }

        # Now safe to prepare the named gene character vector. Here we're
        # passing in the `ensgene` as the value and `symbol` as the name. This
        # works with the constructor function to match the counts matrix by
        # the ensgene, then use the symbol as the name.
        match <- match(x = genes, table = gene2symbol[[format]])
        gene2symbol <- gene2symbol[match, , drop = FALSE]
        ensgene <- gene2symbol[["ensgene"]]
        names(ensgene) <- gene2symbol[["symbol"]]

        .plotGene(
            object = counts,
            genes = ensgene,
            metadata = metadata,
            interestingGroups = interestingGroups,
            color = color,
            countsAxisLabel = countsAxisLabel,
            headerLevel = headerLevel,
            return = return)
    })



#' @rdname plotGene
#' @export
setMethod(
    "plotGene",
    signature("matrix"),
    .plotGene)
